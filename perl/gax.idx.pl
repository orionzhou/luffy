#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gax.idx.pl - index a Gax file (*.bb, *.gz, *.gz.tbi)

=head1 SYNOPSIS
  
  idxgax.pl [-help] [-i input-file] [-s genome-size-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gax)
    -s (--size)   chrom-size file for target genome

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
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;

my ($fi, $fs) = ('', '');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "size|s=s" => \$fs,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fs;

runCmd("sort -k1,1 -k2,2n -k3,3n $fi -o $fi");
runCmd("bgzip -c $fi > $fi.gz");
runCmd("tabix -s 1 -b 2 -e 3 $fi.gz");

runCmd("gax2bed.pl -i $fi -o $fi.bed");
runCmd("bedSort $fi.bed $fi.bed");
runCmd("bedToBigBed -tab $fi.bed $fs $fi.bb");

