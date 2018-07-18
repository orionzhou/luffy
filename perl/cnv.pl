#!/usr/bin/perl
use strict;
use FindBin;
use lib "$FindBin::Bin";
use InitPath;
use Common;
use Bio::Seq;
use WindowStat;
use Time::HiRes qw/gettimeofday tv_interval/;
use Data::Dumper;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
my $refDb = "mt_35";
my $chr = "chr5";

my $dir = "$ENV{'misc2'}/cnv";
my $opt = "acc26";
my $ids = get_acc_ids($opt);

my $f01  = file($dir, "01_windows.tbl");
my $f02  = file($dir, "02_cov.tbl");

write_cov_windows($f01, $f02, $refDb, $chr, $opt);
#runCmd("cov_window -i $f01 -o $f02 -t $opt -c 2", 1);
sub write_cov_windows {
  my ($f01, $f02, $refDb, $chr, $opt) = @_;
  my $fh = new IO::File $f01, "w";
  print $fh join("\t", qw/id chr beg end/)."\n";
  my ($winSize, $winStep) = (2000, 2000);
  my $locAry = getWindows(-chr=>$chr, -winsize=>$winSize, -winstep=>$winStep, -db=>$refDb, -opt=>3);
  for my $i (0..@$locAry-1) {
    my ($wbeg, $wend) = @{$locAry->[$i]};
    print $fh join("\t", $i+1, $chr, $wbeg, $wend)."\n";
  }
  close $fh;
}


