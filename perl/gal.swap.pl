#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal.swap.pl - swap QRY and TGT fields

=head1 SYNOPSIS
  
  galswap.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Location;
use Gal;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('', '');
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

print $fho join("\t", @HEAD_GAL)."\n";
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps == 20;
  my ($id, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $ali, $mat, $mis, $qN, $tN, $ident, $score, $tLocS, $qLocS) = @$ps;
  $tSrd eq "+" || die "$id: $tId tSrd: $tSrd\n";
  my ($qloc, $tloc) = (locStr2Ary($qLocS), locStr2Ary($tLocS));
  @$qloc == @$tloc || die "unequal pieces\n";

  if($qSrd eq "+") {
    ($tLocS, $qLocS) = ($qLocS, $tLocS);
  } else {
    my $tlen = $tEnd - $tBeg + 1;
    my $qlen = $qEnd - $qBeg + 1;
    $tloc = [ reverse map {[$tlen-$_->[1]+1, $tlen-$_->[0]+1]} @$tloc ];
    $qloc = [ reverse map {[$qlen-$_->[1]+1, $qlen-$_->[0]+1]} @$qloc ];
    ($tLocS, $qLocS) = (locAry2Str($qloc), locAry2Str($tloc));
  }
  
  ($tId, $qId) = ($qId, $tId);
  ($tBeg, $qBeg) = ($qBeg, $tBeg);
  ($tEnd, $qEnd) = ($qEnd, $tEnd);
  ($tSize, $qSize) = ($qSize, $tSize);
  ($tN, $qN) = ($qN, $tN);
  
  print $fho join("\t", $id, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $ali, $mat, $mis, $tN, $qN, $ident, $score, $tLocS, $qLocS)."\n";
}
close $fhi;
close $fho;

__END__
