#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  prb.expand.pl - expand probe location

=head1 SYNOPSIS
  
  prb.expand.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (PRB)
    -o (--out)    output file (BED)
    -d (--dist)   probe distance (default: 150)
    -n (--one)    only expand to one direction (default: off)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/min max sum/;

#--------------------------------- MAIN -----------------------------------#
my ($fi, $fo) = ('', '');
my $dist = 150;
my $one_flag;
my $help_flag;
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "dist|d=i" => \$dist,
  "one|n"  => \$one_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my ($fhi, $fho);
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

while(<$fhi>) {
  chomp;
  my ($chr, $pos) = split "\t";
  my $beg = $one_flag ? $pos : max(1, $pos - $dist);
  my $end = $pos + $dist;
  print $fho join("\t", $chr, $beg - 1, $end)."\n";
}
close $fhi;
close $fho;

exit 0;
