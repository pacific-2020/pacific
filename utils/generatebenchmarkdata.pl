#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
my $dataid = shift;
my %coverage = (
	"Cornidoviridae.genomes.fasta" => 47.5,
	"Homo_sapiens.GRCh38.cdna.e99.canonicaltranscripts.fasta" => 0.227,
	"Influenza.genomes.fasta" => 9.65,
	"Metapneumovirus.genomes.fasta" => 228,
	"Rhinovirus.genomes.fasta" => 16.7,
	"Sars_cov_2.genomes.fasta" => 5.8,
	"Unrelated.genomes.fa" => 0.01738
);


for (my $i = $dataid; $i <= $dataid; $i++){
	my $seed = $dataid + 2020;
	print "generating $i test\n";
	open (OUT, ">testdata$i.fasta");
	foreach my $class (sort keys %coverage){
		system("~/art_bin_MountRainier/art_illumina -rs $seed -ss HS25 -sam -M -o test$i.$class -l 150 -f $coverage{$class} -i $class");
		my %reads = ();
		open (F, "<test$i.$class.sam") or die $!;
		while (<F>){
			chomp $_;
			next if (substr($_,0,1) eq "@");
			my @a = split ("\t", $_);
			$reads{$a[0]} = "$a[3]:$a[1]:" . (($a[5] ne "150M") ? "INDEL" : "N");
		}
		close F;
		open (F, "<test$i.$class.aln") or die $!;
		while (<F>){
			chomp $_;
			if ($_ =~ />\S+\s+(\S+)/){
				my $rid = $1;
				my $refseq = <F>;
				my $readseq = <F>;
				if ($refseq ne $readseq ){
					$readseq =~ s/-//g;
					if ($reads{$rid} =~ /:N$/){
						$reads{$rid} =~ s/:N$/:MM/;
						print OUT ">T$i:$rid:$reads{$rid}\n$readseq";
					}
					else {
						print OUT ">T$i:$rid:$reads{$rid}\n$readseq";
					}
				}
				else {
					$reads{$rid} =~ s/:N$/:EQUAL/;
					print OUT ">T$i:$rid:$reads{$rid}\n$readseq";
				}
			}
		}
		close F;		
	}
	close OUT;
}
exit;
