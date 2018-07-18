#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  bedfilter.pl - filter a BED file

=head1 SYNOPSIS
  
  bedfilter.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (BED format)
    -o (--out)    output file (default: stdout)
    -l (--len)    minimum length (default: 1) 

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;

my ($fi, $fo) = ('') x 2;
my $minlen = 1;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "len|l=i" => \$minlen,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my ($fhi, $fho);
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $sum = 0;
while(<$fhi>) {
  chomp;
  next if /^\s*$/;
  my ($id, $beg, $end) = split "\t";
  $end - $beg >= $minlen || next;
  print $fho join("\t", $id, $beg, $end)."\n";
}
close $fhi;
close $fho;

exit 0;
