#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  winslide.pl - create sliding windows for given intervals

=head1 SYNOPSIS
  
  winslide.pl [-help] [-in input] [-step win-step] [-size win-size] [-out output]

  Options:
    -help   brief help message
    -in     input file
    -out    output file
    -step   sliding window step (default: 5)
    -size   sliding window size (default: 60)
    -opt    option(1 - strict, 2 - loose; default: 1)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use POSIX qw/ceil floor/;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my $opt = 1;
my ($size, $step) = (60, 5);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "step|t=i" => \$step,
  "size|z=i" => \$size,
  "opt|p=i"  => \$opt
) or pod2usage(2);
pod2usage(1) if $help_flag;

my ($fhi, $fho);
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

#print $fho join("\t", qw/chr beg end/)."\n";
while(<$fhi>) {
  chomp;
  next if /(^#)|(^id\s)|(^chr\s)/;
  my @ps = split "\t";
  my ($chr, $beg, $end);
  next unless @ps >= 2;
  if(@ps >= 3) {
    ($chr, $beg, $end) = @ps;
  } else {
    ($chr, $beg, $end) = ($ps[0], 1, $ps[1]);
  }
  my $wins = sliding_windows($beg, $end, $step, $size);
  for (@$wins) {
    print $fho join("\t", $chr, $_->[0], $_->[1])."\n";
  }
}
close $fhi;
close $fho;

sub sliding_windows {
  my ($beg, $end, $step, $size, $opt) = @_;
  my $n_win = int(($end-$beg+1-$size)/$step) + 1;
  return [] if $n_win < 1;

  my @wins;
  for my $i (0..$n_win-1) {
    my $beg1 = $beg + $step * $i;
    my $end1 = $beg + $step * $i + $size - 1;
    $end1 = min($end1, $end);
    $end1 <= $end || die "[$beg-$end] range error: [$beg1-$end1]\n";
    push @wins, [$beg1-1, $end1];
  }
  return \@wins;
}

exit 0;
