#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal.tiling.pl - tiles a Gal file

=head1 SYNOPSIS
  
  galtiling.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (Gal) file
    -o (--out)    output (Gal) file
    -m (--min)    minimum tiling length (default: 10)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gal;
use Location;
use Data::Dumper;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my $minlen = 10;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "min|m=i" => \$minlen,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

runCmd("tail -n +2 $fi | sort -k7,7 -k8,8n -o $fi.tmp");

my ($fhi, $fho);
open ($fhi, "<$fi.tmp") || die "cannot read $fi.tmp\n";
open ($fho, ">$fo") || die "cannot write $fo\n";
print $fho join("\t", @HEAD_GAL)."\n";

my @lines;
my $cqid = '';
while(<$fhi>) {
  chomp;
  my ($id, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $lev, $ali, $mat, $mis, $qN, $tN, $ident, $score, $tls, $qls) = split "\t";
  my ($tl, $ql) = (locStr2Ary($tls), locStr2Ary($qls));
  my $line = [$id, $tId, $tBeg, $tEnd, $tSrd, $tSize,
      $qId, $qBeg, $qEnd, $qSrd, $qSize, $score, $tl, $ql];
  if($qId eq $cqid) {
    push @lines, $line;
  } elsif($cqid ne '') {
    gal_tile_one(\@lines, $minlen, $fho);
    backOneLine($fhi);
#    printf "%s: %d\n", $cqid, scalar(@lines);
    @lines = ();    
    $cqid = $qId;
  } else {
    push @lines, $line; 
    $cqid = $qId;
  }
}
gal_tile_one(\@lines, $minlen, $fho) if @lines > 0;

close $fhi;
close $fho;
runCmd("rm $fi.tmp");


sub gal_tile_one {
  my ($lines, $minlen, $fho) = @_;
  my @locs = map {[$_->[7], $_->[8]]} @$lines;
  my @scores = map {$_->[-3]} @$lines;
  my $ref = tiling(\@locs, \@scores, 2);
  for (@$ref) {
    my ($qb, $qe, $idx, $idxs) = @$_;
    next if $qe - $qb + 1 < $minlen;
    my ($id, $tId, $tBeg, $tEnd, $tSrd, $tSize, $qId, $qBeg, $qEnd, 
      $qSrd, $qSize, $score, $tl, $ql) = @{$lines->[$idx]};
    my ($rqb, $rqe) = $qSrd eq "+" ? ($qb-$qBeg+1, $qe-$qBeg+1) :
      ($qEnd-$qe+1, $qEnd-$qb+1);

    my ($idx1, $idx2);
    ($idx1, $rqb) = find_interval($rqb, $ql);
    ($idx2, $rqe) = find_interval($rqe, $ql);
    my $rtb = coordTransform($rqb, $ql, "+", $tl, "+");
    my $rte = coordTransform($rqe, $ql, "+", $tl, "+");
    my ($tb, $te) = $tSrd eq "+" ? ($tBeg+$rtb-1, $tBeg+$rte-1) :
      ($tEnd-$rte+1, $tEnd-$rtb+1);
    
    my ($ntl, $ntll) = posOvlp([[$rtb, $rte]], $tl);
    my ($nql, $nqll) = posOvlp([[$rqb, $rqe]], $ql);
    $ntl = [ map {[$_->[0]-$rtb+1, $_->[1]-$rtb+1]} @$ntl ];
    $nql = [ map {[$_->[0]-$rqb+1, $_->[1]-$rqb+1]} @$nql ];
    my ($ntls, $nqls) = (locAry2Str($ntl), locAry2Str($nql));
    my $ali = locAryLen($ntl);
    my $str = "[$qBeg-$qEnd] [$qb-$qe] [$rqb-$rqe] [$nqls]\n".
      "[$tBeg-$tEnd] [$tb-$te] [$rtb-$rte] [$ntls]\n";
    
    print $fho join("\t", $id, $tId, $tb, $te, $tSrd, $tSize,
      $qId, $qb, $qe, $qSrd, $qSize, 
      '', $ali, ('')x5, $score, $ntls, $nqls)."\n";
  }
}
