#!/usr/bin/perl -w
use strict; use Init; use Localdb; use Run;
use Seq; use Gff; use Crp; use Mapping; use Gtb; use Hmm;
use Data::Dumper; use Path::Class;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Math::Round qw/round/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
my $dirW = dir($DIR_Misc2, "crp");

my $d01 = dir($dirW, "01_in");
my $d04 = dir($dirW, "04_hits_picked");
my $d05 = dir($dirW, "05_models");
my $d50 = dir($dirW, "50_aln");

my $f_seq = "/project/youngn/zhoup/Data/genome/mt_35/41_genome.fa";
my $f01_03 = file($d01, "03_con.fa");
my $f05_21 = file($d05, "21.gtb");
my $f05_22 = file($d05, "22_seq.gtb");

my $f50_01 = file($d50, "01_core.fa");
#get_core_seq($f05_21, $f50_01, $f_seq);
#alignCrp(-fi=>$f05_21, -fseq=>$f50_01, -fhmm=>$f01_03, -out=>"$d50/03");

my $f50_11 = file($d50, "11_gene.fa");
gtb2Fas($f05_21, $f50_11, $f_seq);
alignCrp(-fi=>$f05_21, -fseq=>$f50_11, -fhmm=>$f01_03, -out=>"$d50/13");



