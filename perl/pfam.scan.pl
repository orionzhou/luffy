#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  pfam.scan.pl - scan a fasta file with Pfam-A HMMs

=head1 SYNOPSIS
  
  pfam.scan.pl [-help] [-i input-file] [-o output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (fasta format)
    -o (--out)    output file (*.tsv)

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

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my $fm = "$ENV{'data'}/db/pfam/Pfam-A.hmm";
runCmd("hmmscan --cpu $ENV{'nproc'} -o $fo.1.txt $fm $fi");
runCmd("hmmc2htb.pl -i $fo.1.txt -o $fo.2.htb -m $fm -s $fi");
runCmd("htb.qtile.pl -i $fo.2.htb -o $fo.3.htb");
runCmd("htb.filter.pl -i $fo.3.htb -l 10 -e 0.01 -o $fo.4.htb");
runCmd("cut -f2-4,6,7-9,11-13 $fo.4.htb > $fo");
###runCmd("rm $fo.*.*");


