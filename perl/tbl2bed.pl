#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  tbl2bed.pl - convet tab format to bed format

=head1 SYNOPSIS
  
  tbl2bed.pl [-help] [-in input] [-out output]

  Options:
    -help   brief help message
    -in     input (can be 'stdin')
    -out    output (can be 'stdout')

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
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

while(<$fhi>) {
  chomp;
  next if /^\#/;
  my ($id, $val) = split "\t";
  if($id =~ /^(\w+)[\:\-](\d+)\-(\d+)$/) {
    my ($chr, $beg, $end) = ($1, $2, $3);
    print $fho join("\t", $chr, $beg-1, $beg, $val)."\n";
  } else {
    die "unknown line: $id\t$val\n";
  }
}
close $fhi;
close $fho;

exit 0;
