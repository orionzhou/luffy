#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  sort.header.pl - sort a tabular file with header by coordinate

=head1 SYNOPSIS
  
  sort.header.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (Tabular) file
    -o (--out)    output file
    -f (--fmt)    input file format

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Data::Dumper;

my ($fi, $fo) = ('') x 2;
my $fmt = "";
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "fmt|f=s" => \$fmt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fmt;

my ($idxc, $idxb, $idxe);
if($fmt eq "stb" || $fmt eq "tlc") {
  ($idxc, $idxb, $idxe) = (1, 2, 3);
} elsif($fmt eq "gtb") {
  ($idxc, $idxb, $idxe) = (3, 4, 5);
} elsif($fmt eq "gal") {
  ($idxc, $idxb, $idxe) = (2, 3, 4);
} else {
  die "unsupported fmt: $fmt\n";
}

runCmd("(head -n 1 $fi && tail -n +2 $fi | \\
  sort -k$idxc,$idxc -k$idxb,${idxb}n -k$idxe,${idxe}n) > $fo");

exit 0;
