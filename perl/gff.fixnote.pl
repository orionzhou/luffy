#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use Getopt::Long;
use Pod::Usage;
use Common;
use List::MoreUtils qw/first_index/;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

open(FHI, "<$fi") || die "Can't open file $fi for reading";
open(FHO, ">$fo") || die "Can't open file $fo for writing";

print FHO "##gff-version 3\n";
while(<FHI>) {
  chomp;
  next if !$_ || /^\#/;
  my @ps = split "\t";
  my @ary = split(";", $ps[8]);
  my $idx = first_index {/Note/} @ary;
  if($idx > -1 && $idx < $#ary) {
    @ary = @ary[0..$idx];
    $ary[$idx] =~ s/^Note /Note=/;
    $ps[8] = join(";", @ary);
  }
  print FHO join("\t", @ps)."\n";
}
close FHI;
close FHO;


