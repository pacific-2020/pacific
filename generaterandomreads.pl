#!/usr/bin/perl
use strict;
use warnings;

#hardip@rstudio-nci:~/pacific$ # /media/dstore/covid19/seqid.vs.genomeinfo.txt 
#hardip@rstudio-nci:~/pacific$ # /media/dstore/covid19/virus.genbank.20200324.fa
#V1	V2	V3	taxid	assemblyname	organism	lineage
#GCA_000837105.1	ViralMultiSegProj14079	L14461.1	223359	ViralMultiSegProj14079	Tomato mottle virus-[Florida] (viruses)	Viruses; Geminiviridae; Begomovirus; Tomato mottle virus
#GCA_000837105.1	ViralMultiSegProj14079	L14460.1	223359	ViralMultiSegProj14079	Tomato mottle virus-[Florida] (viruses)	Viruses; Geminiviridae; Begomovirus; Tomato mottle virus

my ($infofile, $seqfile, $outdir, $depth, $maxseqlen) = @ARGV;

my %ginfo = ();
open (F, "<$infofile") or dir $!; ##genomeinfo.txt table
while (<F>){
 chomp $_;
 my @a = split ("\t", $_);
 next if ($a[0] eq "V1");
 $ginfo{"$a[3].$a[0]"."_$a[1]"}{$a[2]} = "";
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

foreach my $g (keys %ginfo){
 open (O, "| gzip >$outdir/$g.d$depth.l$maxseqlen.fa.gz") or die $!;
 foreach my $s (sort keys %{$ginfo{$g}}){
  my $maxstart = length($seqinfo{$s}) - $maxseqlen + 1;
  my $maxbp = length($seqinfo{$s}) * $depth;
  my $bpcount = 0;
  while ($bpcount < $maxbp) {
   my $thisstart = int(rand($maxstart));	  
   my $strand = ($thisstart % 2 == 0) ? 1 : -1; 
   my $thisseq = uc(substr($seqinfo{$s}, $thisstart, $maxseqlen));
   if ($strand < 0){
     $thisseq = reverse ($thisseq);
     $thisseq =~ tr/ACGT/TGCA/;
   }
   print O ">$s:$thisstart:$strand\n$thisseq\n";
   $bpcount += $maxseqlen;
  }
 }
 close O;
}

exit;


