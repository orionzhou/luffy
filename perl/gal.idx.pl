#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal.idx.pl - index a Gal file (*.bb, *.gz, *.gz.tbi)

=head1 SYNOPSIS
  
  gal.idx.pl [-help] [-i input-file] [-s genome-size-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gal format)
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

runCmd("gal2psl.pl -i $fi -o $fi.1.psl");
runCmd("pslToBed $fi.1.psl $fi.2.bed");
runCmd("fixbedbygal.pl -i $fi.2.bed -g $fi -o $fi.3.bed");
runCmd("bedSort $fi.3.bed $fi.3.bed");
runCmd("bedToBigBed -tab $fi.3.bed $fs $fi.bb");

runCmd("tail -n +2 $fi | sort -k2,2 -k3,3n -k4,4n -o $fi.4.gal");
runCmd("bgzip -c $fi.4.gal > $fi.gz");
runCmd("tabix -f -s 2 -b 3 -e 4 $fi.gz");

runCmd("rm $fi.1.psl $fi.2.bed $fi.3.bed $fi.4.gal");

__END__

