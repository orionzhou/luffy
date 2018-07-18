#!/usr/bin/perl -w
use strict;
use InitPath;
use Common;
use Path::Class;
use Gff;
use Crp;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
my $refDb = "mt_35";
my $dirW = dir($DIR_Misc2, "crp");


my $f04 = file($dirW, "04_signalp.txt");
#runSigP(-in=>$f03, -out=>$f04, -fedb=>$feDb, -refdb=>$refDb);
my $f05 = file($dirW, "05_groups.txt");
#addSigP(-in1=>$f03, -in2=>$f04, -out=>$f05);
my $f06 = file($dirW, "06_groups.txt");
#addExp(-in=>$f05, -out=>$f06, -ds=>"rnaseq");
my $f07 = file($dirW, "07_groups.txt");
#addEvalue(-in1=>$f06, -in2=>$f01, -out=>$f07);


my $d11 =  dir($dirW, "11_fgenesh_data");
#predictGene1(-in=>$f01, -outdir=>$d11, -refdb=>$refDb);
my $d12 =  dir($dirW, "12_fgenesh_result");
my $f13 = file($dirW, "13_fgenesh.txt");
#predictGene2(-in=>$f01, -dir=>$d12, -out=>$f13, -refdb=>$refDb);
my $f14 = file($dirW, "14_pred.txt");
#predictGene3(-in=>$f13, -out=>$f14);
my $f15 = file($dirW, "15_pred_sum.txt");
#predictGene4(-in1=>$f01, in2=>$f14, -out=>$f15, -refdb=>$refDb);
my $d16 = file($dirW, "16_aln");
#predictGene5(-in=>$f15, -dir=>$d16, -refdb=>$refDb);
#manually examine fgenesh+ predictions
my $f17 = file($dirW, "17_pred_sum.txt");
my $f18 = file($dirW, "18_pred_sum.txt");
#predictGene6(-in=>$f17, -out=>$f18, -refdb=>$refDb);
my $f19 = file($dirW, "19_pred.gff");
my $f20 = file($dirW, "20_extend.gff");
#predictGene7(-in=>$f18, -out1=>$f19, -out2=>$f20);
my $f19b = file($dirW, "19_pred_2.gff");
my $f20b = file($dirW, "20_extend_2.gff");
#convGffLoc(-in=>$f19, -out=>$f19b, -refdb=>$refDb);
#convGffLoc(-in=>$f20, -out=>$f20b, -refdb=>$refDb);

my $d30 =  dir($dirW, "30_gff");
#exportGffs(-in=>$f01, -outdir=>$d30, -opt=>$opt);
my $f32 = file($DIR_Work, '32_hits.gff');
#table2Gff(-in=>$f31, -out=>$f32, -type=>'CDS_supported_by_domain_match', -source=>'hmmsearch', -col_id=>"hitId", -col_loc=>'location', -col_score=>'e_value', -col_note=>['status', 'annotation', 'mt_35_id']);


