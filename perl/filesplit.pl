#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  filesplit.pl - split a file into N pieces with equal records

=head1 SYNOPSIS
  
  filesplit.pl [-help] [-in input-file] [-out output-prefix]
                       [-n number-files] [-l lines-per-record]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output prefix
    -n (--num)    number of output files (def: 1)
    -l (--line)   lines per record (def: 1)

=head1 BUGS
  
=head1 REFERENCES
  
=head1 VERSION
  
  0.1
  
=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use POSIX qw/ceil floor/;

my ($fi, $fo) = ('') x 2;
my ($n, $l) = (1, 1);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "n=i"      => \$n,
  "l=i"      => \$l,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;
pod2usage(2) if !$n || !$l;

my $lines = runCmd("wc -l $fi", 2);
my @ps = split(" ", $lines->[0]);
my $n_line = $ps[0];
my $n_rec = $n_line / $l;

#$lines = runCmd("grep -c '>' $fi", 2);
#$n_line == $n_rec * 2 || die "n_line[$n_line] != 2 * n_fas[$n_fas]\n";

my $rec_per_file = ceil($n_rec / $n);
my $line_per_file = $rec_per_file * $l;

printf "%d lines in total => %d files [%d lines]\n", $n_line, $n, $line_per_file;
runCmd("split -l $line_per_file -d -a 2 $fi $fo");

exit 0;



