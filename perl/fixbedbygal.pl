#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  fixbedbygal.pl - fix 'name' field in Bed by Gal info

=head1 SYNOPSIS
  
  fixbedbygal.pl [-help] [-i input-Bed] [-g Gal-file] [-o output-Bed]

  Options:
      -h (--help)   brief help message
      -i (--in)     input file (Bed format)
      -g (--gal)    Gal file
      -o (--out)    output (Bed) file

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

my ($fi, $fg, $fo) = ('', '', '');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "gal|g=s"  => \$fg,
  "out|o=s"  => \$fo
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fg || !$fo;

my ($fhi, $fho, $fhg);
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

open ($fhg, "<$fg") || die "Can't open file $fg for writing: $!\n";

my $line = <$fhg>;
while(<$fhi>) {
  chomp;
  my @ps = split "\t";
  my ($chr, $beg, $end) = @ps[0..2];

  $line = <$fhg>;
  my ($id, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize) = split("\t", $line);
  die "sync error $chr:$beg-$end\n" if $chr ne $tId || $beg ne $tBeg-1 || $end ne $tEnd;

  $ps[3] = sprintf("%s[%d]:%d-%d", $qId, $qSize, $qBeg, $qEnd);
  print $fho join("\t", @ps)."\n";
}
close $fhi;
close $fho;
close $fhg;


