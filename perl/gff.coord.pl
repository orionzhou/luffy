#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gff.coord.pl - transform coordinates in the input GFF using existing path file 

=head1 SYNOPSIS
  
  gff.coord.pl [-help] [-in input-file] [-path path-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -p (--path)   path file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Location;

my ($fi, $fo, $fp) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "path|p=s" => \$fp,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fp;

open(FHI, "<$fi") || die "Can't open file $fi for reading";
open(FHO, ">$fo") || die "Can't open file $fo for writing";

my $t = readTable(-in=>$fp, -header=>1);
my $h;
for my $i (0..$t->nofRow-1) {
  my ($id, $chr, $beg, $end, $strd) = $t->row($i);
  die "$id has >1 mappings\n" if exists $h->{$id};
  $h->{$id} = [$chr, $beg, $end, $strd];
}

print FHO "##gff-version 3\n";
while(<FHI>) {
  chomp;
  next if !$_ || /^\#/;
  my @ps = split "\t";
  my ($id, $begL, $endL, $srdL) = @ps[0,3,4,6];
  if(exists $h->{$id}) {
    my ($chr, $begI, $endI, $srdI) = @{$h->{$id}};
    my $locI = [ [1, $endI-$begI+1] ];
    my $locO = [ [$begI, $endI] ];

    my $begG = coordTransform($begL, $locI, $srdI, $locO, "+");
    my $endG = coordTransform($endL, $locI, $srdI, $locO, "+");
    ($begG, $endG) = ($endG, $begG) if $begG > $endG;
    my $srdG = $srdI eq "+" ? $srdL : get_revsrd($srdL);

    $ps[0] = $chr;
    $ps[3] = $begG;
    $ps[4] = $endG;
    $ps[6] = $srdG;
  }
  print FHO join("\t", @ps)."\n";
}
close FHI;
close FHO;
print "you might want to run this:\n";
print "sed -i 's/\\\\//g' $fo\n";

