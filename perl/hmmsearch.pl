#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  hmmsearch.pl - scan a fasta file for HMMs

=head1 SYNOPSIS
  
  hmm.search.pl [-help] [-i input-file] [-m hmm-file] [-o output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (fasta format)
    -o (--out)    output file (Tbl)
    -m (--hmm)    HMM profile

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

my ($fi, $fo, $fm) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "hmm|m=s"  => \$fm,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fm;

runCmd("hmmsearch --domE 10 --cpu 8 -o $fo.1.txt $fm $fi");
runCmd("hmme2htb.pl -i $fo.1.txt -o $fo.2.htb -m $fm -s $fi");
runCmd("htb.filter.pl -i $fo.2.htb -e 0.01 -s 20 -o $fo.3.htb");
runCmd("htb.hmerge.pl -i $fo.3.htb | htb.best.pl | htb.filter.pl -s 50 -q 0.3 -t 0.4 -o $fo.4.htb");
runCmd("cut -f2-4,5,7-9,11-13 $fo.4.htb > $fo");
#runCmd("rm $fo.*.*");
