#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw/shuffle/;

my ($dataid, $testclass, $maxreads, $testproportion) = @ARGV;

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
	shuffleseq($testclass, "tmpi.$i.p$testproportion.$testclass");
	open (OUT, ">fnrtest$i.p$testproportion.$testclass");

	my $classmax = ($maxreads * $testproportion) / 100; ##maximum reads for the proportion
	my $rcount = int($classmax / $seqcounts{$testclass}) + 1; ##per sequence read counts
	system("~/art_bin_MountRainier/art_illumina -q -ss HS25 -sam -M -o tmp.fnrtest$i.p$testproportion.$testclass -l 150 --rcount $rcount -i tmpi.$i.p$testproportion.$testclass");

	my %reads = ();
	open (F, "<tmp.fnrtest$i.p$testproportion.$testclass.sam") or die $!;
	while (<F>){
		chomp $_;
		next if (substr($_,0,1) eq "@");
		my @a = split ("\t", $_);
		$reads{$a[0]} = "$a[3]:$a[1]:" . (($a[5] ne "150M") ? "INDEL" : "N");
	}
	close F;
	my $total = 0;	
	open (F, "<tmp.fnrtest$i.p$testproportion.$testclass.aln") or die $!;
	while (<F>){
		chomp $_;
		if ($_ =~ />\S+\s+(\S+)/){
			my $rid = $1;
			my $refseq = <F>;
			my $readseq = <F>;
			last if ($total >= $classmax);
			$total++;
			if ($refseq ne $readseq ){
				$readseq =~ s/-//g;
				if ($reads{$rid} =~ /:N$/){
					$reads{$rid} =~ s/:N$/:MM/;
					print OUT ">FNRT$i:P$testproportion:$rid:$reads{$rid}\n$readseq";
				}
				else {
					print OUT ">FNRT$i:P$testproportion:$rid:$reads{$rid}\n$readseq";
				}
			}
			else {
				$reads{$rid} =~ s/:N$/:EQUAL/;
				print OUT ">FNRT$i:P$testproportion:$rid:$reads{$rid}\n$readseq";
			}
		}
	}
	close F;


	shuffleseq("Homo_sapiens.GRCh38.cdna.e99.canonicaltranscripts.fasta", "tmpi.$i.p$testproportion.Homo_sapiens.GRCh38.cdna.e99.canonicaltranscripts.fasta");
	$classmax = $maxreads - $total; ##maximum reads of human
	$rcount = int($classmax / $seqcounts{"Homo_sapiens.GRCh38.cdna.e99.canonicaltranscripts.fasta"}) + 1; ##per sequence read counts
	system("~/art_bin_MountRainier/art_illumina -q -ss HS25 -sam -M -o tmp.fnrtest$i.p$testproportion.$testclass -l 150 --rcount $rcount -i tmpi.$i.p$testproportion.Homo_sapiens.GRCh38.cdna.e99.canonicaltranscripts.fasta");
	%reads = ();
	open (F, "<tmp.fnrtest$i.p$testproportion.$testclass.sam") or die $!;
	while (<F>){
		chomp $_;
		next if (substr($_,0,1) eq "@");
		my @a = split ("\t", $_);
		$reads{$a[0]} = "$a[3]:$a[1]:" . (($a[5] ne "150M") ? "INDEL" : "N");
	}
	close F;
	open (F, "<tmp.fnrtest$i.p$testproportion.$testclass.aln") or die $!;
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
					print OUT ">FNRT$i:$rid:$reads{$rid}\n$readseq";
				}
				else {
					print OUT ">FNRT$i:$rid:$reads{$rid}\n$readseq";
				}
			}
			else {
				$reads{$rid} =~ s/:N$/:EQUAL/;
				print OUT ">FNRT$i:$rid:$reads{$rid}\n$readseq";
			}
		}
	}
	close F;
	close OUT;
	unlink("tmp.fnrtest$i.p$testproportion.$testclass.aln");
	unlink("tmp.fnrtest$i.p$testproportion.$testclass.fq");
	unlink("tmp.fnrtest$i.p$testproportion.$testclass.sam");
	unlink("tmpi.$i.p$testproportion.$testclass");
	unlink("tmpi.$i.p$testproportion.Homo_sapiens.GRCh38.cdna.e99.canonicaltranscripts.fasta")
}
exit;


sub shuffleseq {
	##shuffle order of sequences in the file
	my $inputfile = shift;
	my $outputfile = shift;
	my @headers = ();
	my %sequences = ();
	my $seqid = "";
	open (F, "<$inputfile") or die $!;
	while (<F>){
		chomp $_;
		if ($_=~/(>\S+)/){
			$seqid = $1;
			push (@headers, $seqid);
		}
		else {
			$sequences{$seqid} .= $_;
		}
	}
	close F;
	@headers = shuffle @headers;
	open (O, ">$outputfile") or die $!;
	map { print O "$_\n$sequences{$_}\n" } @headers;
	close O;
}

