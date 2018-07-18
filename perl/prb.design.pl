#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  prb.design.pl - design probes for given regions

=head1 SYNOPSIS
  
  prb.design.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (BED)
    -o (--out)    output file (PRB)
    -d (--dist)   minimum probe distance (default: 500)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;

#--------------------------------- MAIN -----------------------------------#
my ($fi, $fo) = ('', '');
my $dist = 500;
my $help_flag;
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "dist|d=i" => \$dist,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my ($fhi, $fho);
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my ($prb_chr, $prb_pos) = ("chr0", -1000);
while(<$fhi>) {
  chomp;
  my ($chr, $beg, $end) = split "\t";
  $beg += 1;
  if( $prb_chr ne $chr ) {
    $prb_chr = $chr;
    $prb_pos = -1000;
  }
  while($prb_pos + $dist < $end) {
    my $bound = $prb_pos + $dist;
    my $pos = $bound < $beg ? $beg : $bound + 1;
    $prb_pos = $pos;
    print $fho join("\t", $chr, $pos)."\n";
  }
}

close $fhi;
close $fho;

exit 0;
