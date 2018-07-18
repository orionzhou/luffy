#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  psl2gal.pl - convert a PSL file to GAL format

=head1 SYNOPSIS
  
  psl2gal.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (PSL)
    -o (--out)    output file (Gal)

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

print $fho join("\t", @HEAD_GAL)."\n";

my $id = 1;
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
  my $alnLen = $match + $misMatch + $repMatch + $baseN;
  $alnLen == sum(@blockSizes) || die "block size error:$alnLen/".sum(@blockSizes)."\n";
  $alnLen + $qBaseIns == $qEnd-$qBeg || die "qLen error\n$qId $tId\n";
  $alnLen + $tBaseIns == $tEnd-$tBeg || die "hLen error\n";
  
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
  my ($mat, $mis, $qN, $tN) = ($match, $misMatch, $baseN, 0);
  my $score = $mat - $mis - $qNumIns - $tNumIns;
  my $ali = $mat + $mis + $qN + $tN;
  my $ident = sprintf "%.03f", $mat / ($mat + $mis);
  print $fho join("\t", $id, $tId, $tBeg+1, $tEnd, $tSrd, $tSize,
    $qId, $qBeg+1, $qEnd, $qSrd, $qSize, 
    '', $ali, $mat, $mis, $qN, $tN, $ident, $score, $tLocS, $qLocS)."\n";
  $id += 1;
}
close $fhi;
close $fho;


__END__
