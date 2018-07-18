#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use Getopt::Long;
use Pod::Usage;
use Common;

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
  my $h = parse_gff_tags($ps[8]);
  
  if(exists $h->{"conf_class"}) {
    $h->{"Note"} = sprintf("[%s]%s", $h->{"conf_class"}, $h->{"Note"});
    $ps[8] = join(";", map {$_."=".$h->{$_}} keys(%$h));
  }
 
#  my $tag_te = $ps[2] eq "gene" && exists $h->{"Name"} && $h->{"Name"} =~ /Medtr\w{1,2}te/;
#  $ps[2] = "transposable_element_gene" if $tag_te;
  
  $ps[2] = "gene" if $ps[2] eq "transposable_element";
#    $ps[7] = ".";
  print FHO join("\t", @ps)."\n";
}
close FHI;
close FHO;


