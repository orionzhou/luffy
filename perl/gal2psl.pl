#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal2psl.pl - convert GAL file to PSL format

=head1 SYNOPSIS
  
  gal2psl.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gal)
    -o (--out)    output file (PSL)

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

while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps == 21;
  my ($cid, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $lev, $ali, $mat, $mis, $qN, $tN, $ident, $score, $tlS, $qlS) = @$ps;
  my ($qLoc, $tLoc) = (locStr2Ary($qlS), locStr2Ary($tlS));
  $tSrd eq "+" || die "$cid: tSrd -\n";

  @$qLoc == @$tLoc || die "unequal pieces\n";
  my $nBlock = @$qLoc;
  
  my (@blockSizes, @qBegs, @tBegs);
  my (@qIns, @tIns);
  my ($rqe_p, $rte_p);
  for my $i (0..$nBlock-1) {
    my ($rqb, $rqe) = @{$qLoc->[$i]};
    my ($rtb, $rte) = @{$tLoc->[$i]};
    my ($len, $len2) = ($rqe-$rqb+1, $rte-$rtb+1);
    $len == $len2 || 
      die "block size unequal: $qId-$tId $rqb-$rqe : $rtb-$rte\n";
    my $tb = $tBeg + $rtb - 1;
    my $qb = $qSrd eq "-" ? $qSize-$qEnd+1 + $rqb-1 : $qBeg + $rqb - 1;
    
    push @blockSizes, $len;
    push @tBegs, $tb-1;
    push @qBegs, $qb-1;
    if($i > 0) {
      my $tIns = $rtb - $rte_p - 1;
      my $qIns = $rqb - $rqe_p - 1;
      push @tIns, $tIns if $tIns > 0;
      push @qIns, $qIns if $qIns > 0;
    }
    ($rqe_p, $rte_p) = ($rqe, $rte);
  }
  my ($qNumIns, $tNumIns) = (scalar(@qIns), scalar(@tIns));
  my ($qBaseIns, $tBaseIns) = (0, 0);
  $qBaseIns = sum(@qIns) if $qNumIns > 0;
  $tBaseIns = sum(@tIns) if $tNumIns > 0;
  my $blockSizes = join(",", @blockSizes).",";
  my $qBegs = join(",", @qBegs).",";
  my $tBegs = join(",", @tBegs).",";
  
  my $repMatch = 0;
  ($mat, $mis) = ($ali, 0) if $mat eq "";
  $qN ||= 0;
  $tN ||= 0;
  print $fho join("\t", $mat, $mis, $repMatch, $qN+$tN, 
    $qNumIns, $qBaseIns, $tNumIns, $tBaseIns, $qSrd, 
    $qId, $qSize, $qBeg-1, $qEnd, $tId, $tSize, $tBeg-1, $tEnd, 
    $nBlock, $blockSizes, $qBegs, $tBegs)."\n";
}
close $fhi;
close $fho;

__END__
