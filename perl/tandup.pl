#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use Common;
use Time::HiRes qw/gettimeofday tv_interval/;

my $dirW = dir($DIR_Misc2, "tandup");
my $d01 =  dir($dirW, "01_gene_family");
my $f02 = file($dirW, "02_gene_family.txt");
#mergeClusterFiles(-dir=>$d01, -out=>$f02, -opt=>1);
my $d11 =  dir($dirW, "11_gene_family_lc");
my $f12 = file($dirW, "12_gene_family_lc.txt");
#mergeClusterFiles(-dir=>$d11, -out=>$f12, -opt=>2);
my $f21 = file($dirW, "21_lc.txt");

sub mergeClusterFiles {
  my ($dir, $fo, $opt) = rearrange(['dir', 'out', 'opt'], @_);
  opendir(DIRI, $dir) or die "cannot open $dir\n";
  $opt ||= 1;
  my $header = [qw/cid gid chr gcnt/];
  $header = [qw/cid lcid gid chr gcnt/] if $opt == 2;
  my $t = Data::Table->new([], $header, 0);
  for my $fn (sort readdir(DIRI)) {
    next unless $fn =~ /^cl(\d+)/;
    my $cid = $1;
    my $fi = file($dir, $fn);
    my $fh = new IO::File $fi, "r";
    my $lcid = 0;
    while(<$fh>) {
      chomp;
      next unless $_;
      my @ps = split "\t";
      if($opt == 1) {
        $t->addRow([$cid, @ps[2, 0, 1]]);
      } elsif($opt == 2) {
        $lcid ++ if @ps == 3;
        $t->addRow([$cid, $lcid, @ps[2, 0, 1]]);
      }
    }
  }
  $t->sort("cid", 0, 0, "gid", 1, 0) if $opt == 1;
  $t->sort("cid", 0, 0, "lcid", 0, 0, "gid", 1, 0) if $opt == 2;
  my $fh = new IO::File $fo, "w";
  print $fh $t->tsv(1);
}


