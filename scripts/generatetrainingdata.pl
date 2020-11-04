#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my ($inputfasta, $maxseqlen, $depth, $strand, $outputfasta) = @ARGV;
GetOptions ("inputfasta|i=s" => \$inputfasta,    # path to input fasta file for sampling
            "maxseqlen|l=i"   => \$maxseqlen,      # numeric integer for maximum length of the sequence
            "depth|d=i"  => \$depth,
            "strand|s=s"  => \$strand,
         		"outputfasta|o=s" => \$outputfasta,    # path to output fasta file 
            );


my $seqinfo = readFasta($inputfasta);

if ($outputfasta =~ /gz$/){
	open (O, "| gzip >$outputfasta") or die $!;
}
else{
	open (O, ">$outputfasta") or die $!;
}
	
foreach my $s (sort keys %$seqinfo) {
	my $maxstart = length($$seqinfo{$s}) - $maxseqlen + 1;
	my $maxbp = length($$seqinfo{$s}) * $depth;
	my $bpcount = 0;
	while ($bpcount < $maxbp) {
		my $thisstart = int(rand($maxstart));	  
		my $thisseq = uc(substr($$seqinfo{$s}, $thisstart, $maxseqlen));
		my $thisstrand = 1;
		if ($strand eq "both"){
			$thisstrand = ($thisstart % 2 == 0) ? 1 : -1; 
			if ($thisstrand < 0){
				$thisseq = reverse ($thisseq);
				$thisseq =~ tr/ACGT/TGCA/;
			}
		}
		elsif ($strand eq "antisense"){
			$thisseq = reverse ($thisseq);
			$thisseq =~ tr/ACGT/TGCA/;
		}	
		print O ">$s:$thisstart:$thisstrand\n$thisseq\n";
		$bpcount += $maxseqlen;
	}
}
close O;

sub readFasta {
  my $file = shift;
  my %s = ();
  my $header = "";
  if ($file =~ /gz$/){
  	open (F, "gunzip -c $file |") or die $!; 
  }
  else{
  	open (F, "<$file") or die $!;
  }
  while (<F>){
    chomp $_;
    if ($_ =~ />(\S+)/){
      $header = $1;
    }
    else {
      $s{$header} .= $_;
    }
  }
  close F;
  return \%s;
}

