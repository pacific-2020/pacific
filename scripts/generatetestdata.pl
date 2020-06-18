#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/shuffle/;

#hardip@rstudio-nci:~/pacific$ # /media/dstore/covid19/seqid.vs.genomeinfo.txt 
#hardip@rstudio-nci:~/pacific$ # /media/dstore/covid19/virus.genbank.20200324.fa
#V1	V2	V3	taxid	assemblyname	organism	lineage
#GCA_000837105.1	ViralMultiSegProj14079	L14461.1	223359	ViralMultiSegProj14079	Tomato mottle virus-[Florida] (viruses)	Viruses; Geminiviridae; Begomovirus; Tomato mottle virus
#GCA_000837105.1	ViralMultiSegProj14079	L14460.1	223359	ViralMultiSegProj14079	Tomato mottle virus-[Florida] (viruses)	Viruses; Geminiviridae; Begomovirus; Tomato mottle virus

my ($infofile, $seqfile, $selectionfile, $outdir, $maxreads, $maxseqlen, $totaldset) = @ARGV;

my %ginfo = ();
open (F, "<$infofile") or dir $!; ##genomeinfo.txt table
while (<F>){
	chomp $_;
	my @a = split ("\t", $_);
	next if ($a[0] eq "V1");
	push(@{$ginfo{"$a[3].$a[0]"."_$a[1]"}}, $a[2]);
}
close F;
print "done reading info file\n";

my %seqinfo = ();
my $sid = "";
open (F, "<$seqfile") or die $!; ##genomic sequences
while (<F>) {
	chomp $_;
	if ($_ =~ />(\S+)/) {
		$sid = $1;
	}
	else{
		$seqinfo{$sid} .= $_;
	}
}
close F;

print "done reading fasta file\n";

my %selections = ();
my @groups = ();
open (F, "<$selectionfile") or die $!; ##a file of the format Cornidovirineae:28295.GCA_900197455.1_PEDV_GER_L01061-K07_15-03_2015.d30.l150.fa
while (<F>){
	chomp $_;
	$_=~/(\S+):(\d+\.G\S+).d30.l150.fa/;
	my $group = $1;
	my $genome = $2;
	push(@{$selections{$group}},  $genome);
}
close F;
##add selections for random
@{$selections{'random'}} = ('random');

map { push (@groups, $_) } keys %selections;

my @characters = ("A", "C", "G", "T");

open (X, ">$outdir/tdi.info.txt") or die $!;
for (my $id = 1; $id<= $totaldset; $id++){
	print "Generating tid $id data\n"; 
	open (O, "| gzip >$outdir/tdi$id.fq.gz") or die $!;
	my %counter = ();
	my $thisseq = "";
	my $reads2select = 0;
	my $remaining2select = $maxreads;
	my $readsselected = 0;
	@groups = shuffle @groups;
	foreach my $sel (@groups){
		$reads2select = (keys (%counter) < keys (%selections) - 1) ? int(rand($remaining2select)) : $remaining2select;
		$readsselected+=$reads2select;
		$remaining2select = $maxreads - $readsselected;
		#print "#tdi$id,$sel,$reads2select,$remaining2select,".(keys %counter)."\n";
		print X "#tdi$id,$sel,$reads2select\n";
		next if $reads2select == 0; ##move on if there are no reads to select for this group
		if ($sel eq "random"){
			for (my $r=1;$r<=$reads2select;$r++){
				$thisseq = "";
				while (length($thisseq) < $maxseqlen) {
					$thisseq .= $characters[int(rand(4))];
				}
				print O "\@tdi$id:$sel:$sel:$r/$reads2select\n$thisseq\n+\n".("I" x $maxseqlen)."\n";
				$counter{$sel}{$sel}++;
			}
		}
		else{
			for (my $r=1;$r<=$reads2select;$r++){
				my $g = ${$selections{$sel}}[int(rand(@{$selections{$sel}}))];
				my $s = ${$ginfo{$g}}[int(rand(@{$ginfo{$g}}))];
				my $maxstart = length($seqinfo{$s}) - $maxseqlen + 1;
				my $thisstart = int(rand($maxstart));
				my $strand = ($thisstart % 2 == 0) ? 1 : -1;
        $thisseq = uc(substr($seqinfo{$s}, $thisstart, $maxseqlen));
        if ($strand < 0){
          $thisseq = reverse ($thisseq);
          $thisseq =~ tr/ACGT/TGCA/;
        }
				print O "\@tdi$id:$sel:$g:$r/$reads2select\n$thisseq\n+\n".("I" x $maxseqlen)."\n";
				$counter{$sel}{$g}++;
			}
		}
	}
	foreach my $c (sort keys %counter){
		map { print X "tdi$id,$c,$_,$counter{$c}{$_}\n"} sort keys %{$counter{$c}};
	}
	close O;
}
close X;
exit;

