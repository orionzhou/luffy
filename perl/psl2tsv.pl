#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  psl2tsv.pl - convert a psl file to tsv format

=head1 SYNOPSIS
  
  psl2tsv.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (.psl)
    -o (--out)    output file (.tsv)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my @colnames = qw/
  qId qBeg qEnd qSrd qSize
  tId tBeg tEnd tSrd tSize
  alnLen match misMatch baseN qNumIns tNumIns qBaseIns tBaseIns ident score
  qLoc tLoc/;
print $fho join("\t", @colnames)."\n";

my ($sMatch, $sMisMatch, $sGapOpen, $sGapExtend) = (2, -3, -5, -2);
while(<$fhi>) {
  chomp;
  my $ps = [ split "\t" ];
  next unless @$ps == 21;
  next if /^(psLayout)|(match)/;
  my ($match, $misMatch, $repMatch, $baseN, 
    $qNumIns, $qBaseIns, $tNumIns, $tBaseIns, $qSrd, 
    $qId, $qSize, $qBeg, $qEnd, $tId, $tSize, $tBeg, $tEnd, 
    $blockNum, $blockSizes, $qBegs, $tBegs) = @$ps;
  
  my $tSrd = "+";
  my @qBegs = split(",", $qBegs);
  my @tBegs = split(",", $tBegs);
  my @blockSizes = split(",", $blockSizes);
  
  $blockNum == @qBegs || die "unequal pieces\n";
  $blockNum == @tBegs || die "unequal pieces\n";
  $blockNum == @blockSizes || die "unequal pieces\n";
  $match += $repMatch;
  my $alnLen = $match + $misMatch + $baseN;
  $alnLen == sum(@blockSizes) || die "$qId block size error:$alnLen/".sum(@blockSizes)."\n";
  $alnLen + $qBaseIns == $qEnd-$qBeg || die "$qId qLen error\n";
  $alnLen + $tBaseIns == $tEnd-$tBeg || die "$qId tLen error\n";
  
  my (@qLoc, @tLoc);
  for my $i (0..$blockNum-1) {
    my $len = $blockSizes[$i];
    my $tb = $tBegs[$i] + 1;
    my $te = $tb + $len - 1;
    
    my ($qb, $qe);
    if($qSrd eq "+") {
      $qb = $qBegs[$i] + 1;
      $qe = $qb + $len - 1;
    } else {
      $qSrd eq "-" || die "unknown strand $qSrd\n";
      $qe = $qSize - $qBegs[$i];
      $qb = $qe - $len + 1;
    }
    
    my $rtb = $tBegs[$i] - $tBeg;
    my $rqb = $qSrd eq "-" ? $qBegs[$i]-($qSize-$qEnd) : $qBegs[$i]-$qBeg;
    push @tLoc, [$rtb+1, $rtb+$len];
    push @qLoc, [$rqb+1, $rqb+$len];
  }
  my ($tLocS, $qLocS) = (locAry2Str(\@tLoc), locAry2Str(\@qLoc));
  
  my (@tIns, @qIns);
  for my $i (1..$blockNum-1) {
    my ($qb, $qe) = @{$qLoc[$i]};
    my ($pqb, $pqe) = @{$qLoc[$i-1]};
    my ($tb, $te) = @{$tLoc[$i]};
    my ($ptb, $pte) = @{$tLoc[$i-1]};
    if ($qb != $pqe + 1) {
      push @qIns, $qb - $pqe - 1;
    }
    if ($tb != $pte + 1) {
      push @tIns, $tb - $pte - 1;
    }
  }
  if (@qIns > 0) {
    sum(@qIns) == $qBaseIns || die "qBaseIns $qBaseIns $qId $tId error\n";
  }
  if (@tIns > 0) {
    sum(@tIns) == $tBaseIns || die "tBaseIns $tBaseIns error\n";
  }
 
  my $score_match = $match * $sMatch;
  my $score_misMatch = $misMatch * $sMisMatch;
  my $numIns = $qNumIns + $tNumIns;
  my $score_indel = 0;
  $score_indel = $sGapOpen + ($numIns - 1) * $sGapExtend if $numIns >= 1;
  my $score = $score_match + $score_misMatch + $score_indel;
  my $ident = sprintf "%.03f", $match / ($match + $misMatch);
  print $fho join("\t", 
    $qId, $qBeg+1, $qEnd, $qSrd, $qSize, 
    $tId, $tBeg+1, $tEnd, $tSrd, $tSize,
    $alnLen, $match, $misMatch, $baseN, 
    $qNumIns, $tNumIns, $qBaseIns, $tBaseIns, $ident, $score,
    $qLocS, $tLocS,
    )."\n";
}
close $fhi;
close $fho;


__END__
