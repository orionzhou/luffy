#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gff.rename.pl - rename GFF Sequence IDs

=head1 SYNOPSIS
  
  gff.rename.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -m (--map)    ID mapping file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;

my ($fi, $fo, $fm) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "map|m=s" => \$fm,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fm;

open(my $fhi, "<$fi") || die "cannot read $fi\n";
open(my $fho, ">$fo") || die "cannot write $fo\n";

my $hm;
open(my $fhm, "<$fm") or die "cannot read $fm\n";
while(<$fhm>) {
  chomp;
  my ($nid, $len, $oid) = split "\t";
  die "$oid has >1 mappings\n" if exists $hm->{$oid};
  $hm->{$oid} = $nid;
}
close $fhm;

print $fho "##gff-version 3\n";
while(<$fhi>) {
  chomp;
  next if !$_ || /^\#/;
  my @ps = split "\t";
  my ($oid, $begL, $endL, $srdL) = @ps[0,3,4,6];
  exists $hm->{$oid} || die "$oid not found\n";
  $ps[0] = $hm->{$oid};
  print $fho join("\t", @ps)."\n";
}
close $fhi;
close $fho;

