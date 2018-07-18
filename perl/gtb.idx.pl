#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.idx.pl - index a Gtb file (*.bb, *.gz, *.gz.tbi)

=head1 SYNOPSIS
  
  gtb.idx.pl [-help] [-i input-file] [-s genome-size-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gtb format)
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

runCmd("(head -n 1 $fi && tail -n +2 $fi | sort -k3,3 -k4,4n -k5,5n) \\
  > $fi.sorted");
runCmd("bgzip -c $fi.sorted > $fi.gz");
runCmd("tabix -S 1 -s 3 -b 4 -e 5 $fi.gz");
runCmd("rm $fi.sorted");

runCmd("gtb2bed.pl -i $fi -o $fi.bed");
runCmd("sort -k1,1 -k2,2n -k3,3n $fi.bed -o $fi.bed");
runCmd("bedToBigBed -tab $fi.bed $fs $fi.bb");
runCmd("rm $fi.bed");


