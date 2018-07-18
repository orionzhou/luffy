#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  samplelines.pl - draw serial/random lines from a file

=head1 SYNOPSIS
  
  samplelines.pl [-help] [options] [-in input-file] [-out output] 

  Options:
    -h, --help       brief help message
    -i, --in         input file
    -o, --out        output
    -n, --sample     number of lines to sample (default: 1)
    -p, --opt        sample option (serial / random; default: serial)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use Time::HiRes qw/gettimeofday tv_interval/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my ($fi, $fo) = ('') x 2;
my ($n_sam, $opt) = (0, 'serial');
my $fhi;
my $fho;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "sample|n=i" => \$n_sam,
  "opt|p=s"    => \$opt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(1) if !$fi;

$fho = \*STDOUT;
unless ($fo eq "stdout" || $fo eq "-" || $fo eq "") {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $n_tot = `wc -l <$fi`;
my @nums;
if($opt eq "serial") {
  @nums = sample_serial($n_tot, $n_sam);
} elsif($opt eq "random") {
  @nums = sample_random($n_tot, $n_sam);
} else {
  die "unknown opt: $opt\n";
}

open ($fhi, "<$fi") || die "cannot read $fi\n";
my ($line, $idx) = (1, 0);
while( <$fhi> ) {
  if($idx < @nums && $line == $nums[$idx]) {
    print $fho $_;
    $idx ++;
  } 
  $line ++;
}
close $fhi;
close $fho;

sub sample_serial {
  my ($n_tot, $n_sam) = @_;
  $n_sam = 1 if $n_sam > $n_tot || $n_sam <= 0;
  my $inc = $n_tot / $n_sam;
  my @nums = map { int($_ * $inc) + 1 } (0..$n_sam-1);
  return @nums;
}
sub sample_random {
  my ($n_tot, $n_sam) = @_;
  $n_sam = 1 if $n_sam > $n_tot || $n_sam <= 0;
  my $inc = $n_tot/$n_sam;
  my @nums = map { int(($_+rand()) * $inc) + 1 } (0..$n_sam-1);
  return @nums;
}

