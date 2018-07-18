#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.rmutr.pl - remove UTRs in a Gtb file

=head1 SYNOPSIS
  
  gtb.rmutr.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gtb)
    -o (--out)    output file (Gtb)

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
use Gtb;

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
if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

print $fho join("\t", @HEAD_GTB)."\n";
while(<$fhi>) {
  chomp;
  /(^id)|(^\#)|(^\s*$)/i && next;
  my $ps = [ split("\t", $_, -1) ];
  @$ps >= 18 || die "not 19 fileds:\n$_\n";
  my ($id, $par, $chr, $beg, $end, $srd, 
    $elocS, $ilocS, $clocS, $flocS, $tlocS, $phase, 
    $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
  if($cat1 ne "mRNA") {
    print "skipped: $id [$cat2]\n";
    next;
  }
  die "$id: no cloc\n" unless $clocS;
  
  my ($rb, $re) = (1, $end - $beg + 1);
  my $cloc = locStr2Ary($clocS);
  $cloc = [ sort {$a->[0] <=> $b->[0]} @$cloc ];
  if($cloc->[0]->[0] > $rb) {
    my $los = $cloc->[0]->[0] - $rb;
    $beg += $los if $srd eq "+";
    $end -= $los if $srd eq "-";
    $cloc = [ map {[$_->[0] - $los, $_->[1] - $los]} @$cloc ];
    ($rb, $re) = (1, $end - $beg + 1);
  }
  if($cloc->[-1]->[1] < $re) {
    my $ros = $re - $cloc->[-1]->[1];
    $end -= $ros if $srd eq "+";
    $beg += $ros if $srd eq "-";
    ($rb, $re) = (1, $end - $beg + 1);
  }
  my ($iloc) = posDiff([[$rb, $re]], $cloc);
 
  $elocS = locAry2Str($cloc);
  $clocS = locAry2Str($cloc);
  $ilocS = locAry2Str($iloc);
  print $fho join("\t", $id, $par, $chr, $beg, $end, $srd, 
    $elocS, $ilocS, $clocS, '', '', $phase, 
    $src, $conf, $cat1, $cat2, $cat3, $note)."\n";
}
close $fhi;
close $fho;



__END__
