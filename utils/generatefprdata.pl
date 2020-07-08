#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw/shuffle/;

my ($dataid, $testclass, $maxreads) = @ARGV;

my %seqcounts = (
	"Coronaviridae.genomes.fasta" => 12,
	"Homo_sapiens.GRCh38.cdna.e99.canonicaltranscripts.fasta" => 22499,
	"Influenza.genomes.fasta" => 1024,
	"Metapneumovirus.genomes.fasta" => 5,
	"Rhinovirus.genomes.fasta" => 130,
	"Sars_cov_2.genomes.fasta" => 87,
	"Unrelated.genomes.fasta" => 65886
);


for (my $i = $dataid; $i <= $dataid; $i++){
	print "generating $i test for $testclass\n";
	open (OUT, ">fprtest$i.$testclass");
	my $total = 0;
	my $currentmax = $maxreads;
	my @files = grep { $_ ne $testclass } keys %seqcounts;
	my $iteration = 0;
	while ($total < $maxreads){
		$iteration++;
		@files = shuffle @files;
		foreach my $class ( @files ){
			my $rcount = int(int(rand($currentmax)) / $seqcounts{$class});
			$rcount = int(($maxreads - $total)/ $seqcounts{$class})  if ($files[$#files] eq $class);
			next if ($rcount == 0 && $iteration == 1);
			$rcount = 1 if ($rcount == 0 && $iteration > 1);
			print "at $iteration, $class, $rcount\n";
			system("~/art_bin_MountRainier/art_illumina -q -ss HS25 -sam -M -o tmp.fprtest$i.$testclass -l 150 --rcount $rcount -i $class");
			my %reads = ();
			open (F, "<tmp.fprtest$i.$testclass.sam") or die $!;
			while (<F>){
				chomp $_;
				next if (substr($_,0,1) eq "@");
				my @a = split ("\t", $_);
				$reads{$a[0]} = "$a[3]:$a[1]:" . (($a[5] ne "150M") ? "INDEL" : "N");
			}
			close F;
			open (F, "<tmp.fprtest$i.$testclass.aln") or die $!;
			while (<F>){
				chomp $_;
				if ($_ =~ />\S+\s+(\S+)/){
					my $rid = $1;
					my $refseq = <F>;
					my $readseq = <F>;
					last if ($total >= $maxreads);
					$total++;
					if ($refseq ne $readseq ){
						$readseq =~ s/-//g;
						if ($reads{$rid} =~ /:N$/){
							$reads{$rid} =~ s/:N$/:MM/;
							print OUT ">FPRT$i:$rid:$reads{$rid}\n$readseq";
						}
						else {
							print OUT ">FPRT$i:$rid:$reads{$rid}\n$readseq";
						}
					}
					else {
						$reads{$rid} =~ s/:N$/:EQUAL/;
						print OUT ">FPRT$i:$rid:$reads{$rid}\n$readseq";
					}
				}
			}
			close F;
			$currentmax = $maxreads - $total;
		}
	}
	close OUT;
	unlink("tmp.fprtest$i.$testclass.aln");
	unlink("tmp.fprtest$i.$testclass.fq");
	unlink("tmp.fprtest$i.$testclass.sam");
}
exit;
