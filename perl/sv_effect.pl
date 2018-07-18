#!/usr/bin/perl -w
use strict; 
use FindBin;
use lib $FindBin::Bin;
use InitPath;
use Path::Class;
use WindowStat;
use Genemodel;

my $refDb = "mt_35";
my $dir = dir($DIR_Repo, $refDb, "40_sv/61_effect");
my $f01 = file($dir, "01.tbl");
my $f02 = file($dir, "02_gene.tbl");
#window_stats_gene(-fi=>$f01, -fo=>$f02, -refdb=>$refDb, -opt=>1);

my $f_gtb = file($DIR_Genome, "mt_35/10_model_Mt3.5v6/62_frame_fixed.gtb");
qry_window_file($f01, $f_gtb, $f02);

