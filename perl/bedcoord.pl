#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  bedcoord.pl - transform the coordinates of a Bed file

=head1 SYNOPSIS
  
  bedcoord.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;

my ($fi, $fo) = ('') x 2;
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

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
  my ($sid, $rb, $re) = split "\t";
  if($sid =~ /^([\w\-\.]+)\-([0-9e\+]+)\-([0-9e\+]+)$/) {
    my ($id, $b, $e) = ($1, $2, $3);
    print $fho join("\t", $id, $b-1 + $rb, $b-1 + $re)."\n";
  } else {
    die "unknown id: $sid\n";
  }
}
close $fhi;
close $fhi;
