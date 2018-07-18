#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  lastz.pl - pipeline of pairwise comparison between 2 genomes

=head1 SYNOPSIS
  
  comp.pl [-help] [-qry query-genome] [-tgt target-genome]

  Options:
    -h (--help)   brief help message
    -q (--qry)    query genome (def: HM056)
    -t (--tgt)    target genome (def: HM101)

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

my ($qry, $tgt) = ('HM056', 'HM101');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "qry|q=s" => \$qry,
  "tgt|t=s" => \$tgt,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $data = $ENV{'data'};
my $qry_fas = "$data/genome/$qry/11_genome.fas";
my $tgt_fas = "$data/genome/$tgt/11_genome.fas";
my $qry_2bit = "$data/db/blat/$qry.2bit";
my $tgt_2bit = "$data/db/blat/$tgt.2bit";
my $qry_size = "$data/genome/$qry/15.sizes";
my $tgt_size = "$data/genome/$tgt/15.sizes";
my $qry_size_bed = "$data/genome/$qry/15.bed";
my $tgt_size_bed = "$data/genome/$tgt/15.bed";
my $qry_gap = "$data/genome/$qry/16_gap.bed";
my $tgt_gap = "$data/genome/$tgt/16_gap.bed";

my $dir = "$data/misc3/$qry\_$tgt/21_lastz";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

run_lastz();
sub run_lastz {
  runCmd("lastz $tgt_fas\[multiple] $qry_fas --format=maf --output=11.maf");
}



__END__

