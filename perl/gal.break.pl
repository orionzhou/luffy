#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal.break.pl - break up long alns into smaller blocks

=head1 SYNOPSIS
  
  gal.break.pl [-help] [-gap gap-size] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -g (--gap)    gap size (default: 1,000)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use Gal;

my ($fi, $fo) = ('') x 2;
my $gap = 1000;
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "gap|g=i"  => \$gap,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

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

print $fho join("\t", @HEAD_GAL)."\n";
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps == 21;
  my ($cid, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $lev, $ali, $mat, $mis, $qN, $tN, $ident, $score, $tlS, $qlS) = @$ps;
  my ($qLoc, $tLoc) = (locStr2Ary($qlS), locStr2Ary($tlS));
  @$qLoc == @$tLoc || die "unequal pieces\n";
  my $nBlock = @$qLoc;
  
  my @idxs = [0, @$qLoc-1];
  for my $i (0..@$qLoc-1) {
    next if $i == 0;
    my ($rqb, $rqe) = @{$qLoc->[$i]};
    my ($rtb, $rte) = @{$tLoc->[$i]};
    my ($prqb, $prqe) = @{$qLoc->[$i-1]};
    my ($prtb, $prte) = @{$tLoc->[$i-1]};
    ($prqe < $rqb && $prte < $rtb) ||
      die "error: $cid $qId\[$prqb-$prqe, $rqb-$rqe] $tId\[$prtb-$prte, $rtb-$rte]\n";
    if($rqb - $prqe - 1 >= $gap) {
      $idxs[-1]->[1] = $i-1;
      push @idxs, [$i, @$qLoc-1];
    }
  }
  
  for my $i (0..@idxs-1) {
    my ($idxb, $idxe) = @{$idxs[$i]};
    my @rql = @$qLoc[$idxb..$idxe];
    my @rtl = @$tLoc[$idxb..$idxe];
    my ($rqb, $rqe) = ($rql[0]->[0], $rql[-1]->[1]);
    my ($rtb, $rte) = ($rtl[0]->[0], $rtl[-1]->[1]);
    my $qb = $qSrd eq "-" ? $qEnd-$rqe+1 : $qBeg+$rqb-1;
    my $qe = $qSrd eq "-" ? $qEnd-$rqb+1 : $qBeg+$rqe-1;
    my $tb = $tSrd eq "-" ? $tEnd-$rte+1 : $tBeg+$rtb-1;
    my $te = $tSrd eq "-" ? $tEnd-$rtb+1 : $tBeg+$rte-1;
    my @nrql = map {[$_->[0]-$rqb+1, $_->[1]-$rqb+1]} @rql;
    my @nrtl = map {[$_->[0]-$rtb+1, $_->[1]-$rtb+1]} @rtl;
    my $qlen = $rqe - $rqb + 1;
    my $tlen = $rte - $rtb + 1;
    my ($qls, $tls) = (locAry2Str(\@nrql), locAry2Str(\@nrtl));
    my $nid = "$cid.".($i+1);
    my $nali = locAryLen(\@nrql);
    print $fho join("\t", $nid, $qId, $qb, $qe, $qSrd, $qSize, 
      $tId, $tb, $te, $tSrd, $tSize, 
      $lev, $nali, $nali, 0, 0, 0, '', '', $qls, $tls)."\n";
  }
}
close $fhi;
close $fho;


__END__
