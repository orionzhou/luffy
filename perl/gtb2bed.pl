#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2bed.pl - convert a Gtb file to BED format

=head1 SYNOPSIS
  
  gtb2bed.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gtb format)
    -o (--out)    output file (Bed format)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Common;
use Location;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my ($fhi, $fho);
if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot read $fi\n";
}
if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $hcol = {
  "default" => "0,0,120",
  "TE"      => "64,64,64",
  "CRP"     => "0,128,255",
  "NBS-LRR" => "255,128,0",
};

#print $fho "#track name=gene_models itemRgb=On useScore=0\n";
while(<$fhi>) {
  chomp;
  /(^id)|(^\#)|(^\s*$)/i && next;
  my $ps = [ split("\t", $_, -1) ];
  @$ps >= 18 || die "not 19 fileds:\n$_\n";
  my ($id, $par, $chr, $beg, $end, $srd, $locES, $locIS, $locCS, $loc5S, $loc3S, $phaseS, $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
  $cat1 eq "mRNA" || next;
  $locCS || die "no CDS for $id\n";
  
  my $idStr = $id;
  $idStr .= "|$note" if $note;
  my $col = exists $hcol->{$cat2} ? $hcol->{$cat2} : $hcol->{"default"};

  my @locs = sort {$a->[0] <=> $b->[0]} @{locStr2Ary($locES)};
  my $n = @locs;
  
  my $rloc = [ sort {$a->[0] <=> $b->[0]} @{locStr2Ary($locCS)} ];
  my ($rtBeg, $rtEnd) = ($rloc->[0]->[0], $rloc->[-1]->[1]);
  my @begs;
  my @lens;
  my ($tBeg, $tEnd);
  if($srd eq "+") {
    @begs = map {$_->[0] - 1} @locs;
    @lens = map {$_->[1] - $_->[0] + 1} @locs;
    ($tBeg, $tEnd) = ($beg+$rtBeg-1, $beg+$rtEnd-1);
  } else {
    @begs = reverse map {$end-$beg+1 - $_->[1]} @locs;
    @lens = reverse map {$_->[1] - $_->[0] + 1} @locs;
    ($tBeg, $tEnd) = ($end-$rtEnd+1, $end-$rtBeg+1);
  }
  print $fho join("\t", $chr, $beg - 1, $end, $idStr, 0, $srd, 
    $tBeg - 1, $tEnd, $col,
    $n, join(",", @lens), join(",", @begs) )."\n";
}
close $fhi;
close $fho;

