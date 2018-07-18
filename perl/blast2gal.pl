#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  blast2gal.pl - convert a BLAST tabular output to GAL format

=head1 SYNOPSIS
  
  blast2gal.pl [-help] [-in input-file] [-out output-file]
    -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen length nident mismatch gaps evalue bitscore qseq sseq'

  Options:
    -help   brief help message
    -in     input file
    -out    output file

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
use Blast;
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
pod2usage(2) if !$fi;

if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

print $fho join("\t", @HEAD_GAL)."\n";

my $id = 1;
while(<$fhi>) {
  chomp;
  my $ps = [ split "\t" ];
  next unless @$ps == 16;
  my ($qId, $qBeg, $qEnd, $qSize, $tId, $tBeg, $tEnd, $tSize, $alnLen, 
    $match, $misMatch, $gaps, $e, $score, $qSeq, $tSeq) = @$ps;
  
  $qBeg <= $qEnd || die "$qId $qBeg > $qEnd\n";
  $alnLen == $match + $misMatch + $gaps 
    || die "len error\n".join("\t", @$ps)."\n";
  my ($qSrd, $tSrd) = ("+") x 2;
  if($tBeg > $tEnd) {
    ($tBeg, $tEnd) = ($tEnd, $tBeg);
    $qSrd = "-";
  }

  my ($qLoc, $tLoc, $stat, $qNumIns, $qIns, $tNumIns, $tIns) = 
    parse_aln_string($qSeq, $tSeq);
  my $mat = sum( map {$_->[0]} @$stat );
  my $mis = sum( map {$_->[1]} @$stat );
  my $qN = sum( map {$_->[2]} @$stat );
  my $tN = 0;
  my $ali = $mat + $mis + $qN + $tN;
  my $nBlock = @$qLoc;
  for my $i (0..$nBlock-1) {
    my ($qbr, $qer) = @{$qLoc->[$i]};
    my ($tbr, $ter) = @{$tLoc->[$i]};
    my ($match, $misMatch, $baseN) = @{$stat->[$i]};
    my $qLen = $qer - $qbr + 1;
    my $tLen = $ter - $tbr + 1;
    $qLen == $tLen || die "len error: $qbr-$qer $tbr-$ter\n";

    my $tb = $qSrd eq "-" ? $tEnd-$ter+1 :$tBeg+$tbr-1;
    my $te = $qSrd eq "-" ? $tEnd-$tbr+1 :$tBeg+$ter-1;
    my $qb = $qBeg + $qbr - 1;
    my $qe = $qBeg + $qer - 1;
  }
  my $ident = $mat / ($mat + $mis);
  my ($tLocS, $qLocS) = (locAry2Str($tLoc), locAry2Str($qLoc));

  print $fho join("\t", $id++, $tId, $tBeg, $tEnd, $tSrd, $tSize,
    $qId, $qBeg, $qEnd, $qSrd, $qSize, 
    '', $ali, $mat, $mis, $qN, $tN, $ident, $score, $tLocS, $qLocS)."\n";
}

__END__
