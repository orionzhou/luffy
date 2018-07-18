#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use File::Path qw/make_path remove_tree/;
use Data::Dumper;
use List::Util qw/min max sum/;
use Common;

my $org = "HM340.APECCA";
my $dir = "/home/youngn/zhoup/Data/genome/$org";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir $dir\n";

runCmd("grep 'NB-ARC' 51.fas.tsv > 42_nbs.tsv");
extract_gtb("40_gene.gtb", "42_nbs.tsv", "42_nbs.gtb");


sub extract_gtb {
  my ($fi, $fl, $fo) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  
  my $hi;
  open(my $fhl, "<$fl") or die "cannot read $fl\n";
  while(<$fhl>) {
    chomp;
    my ($id) = split "\t";
    $hi->{$id} ||= 1;
  }
  printf "%d ids read\n", scalar(keys(%$hi));

  my @idxs;
  for my $i (0..$t->lastRow) {  
    my ($id) = $t->row($i);
    push @idxs, $i if exists $hi->{$id};
  }
  my $ts = $t->subTable(\@idxs);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $ts->tsv(1);
  close $fho;
}
