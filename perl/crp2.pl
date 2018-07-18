#!/usr/bin/perl -w
use strict;
use Common;
use Gtb;
use Crp;
use Genemodel;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $f_seq = "/project/youngn/zhoup/Data/genome/mt_35/41_genome.fa";
my $f_gtb = "/project/youngn/zhoup/Data/genome/mt_35/10_model_Mt3.5v5/62_phase_fixed.gtb";

my $dir = "/project/youngn/zhoup/Data/misc2/crp/05_models";
my $f_hit = "$dir/../04_hits_picked/21_hits.tbl";

my $d01 = "$dir/01_mt35_v5";
my $f01_01 = "$d01/01_eval.tbl";
#compare_models_2(-in=>$f_hit, -gtb=>$f_gtb, -out=>$f01_01);
my $f01_11 = "$d01/11.gtb";
#extract_compatible_models($f01_01, $f_gtb, $f01_11);

my $d05 = "$dir/05_ks";
my $f05_11 = "$d05/11.gtb";

my $d09 = "$dir/09_simple";
my $f09_11 = "$d09/11.gtb";
#build_simple_model($f_hit, $f09_11);

my $f11 = "$dir/11_picked.tbl";
my $f31a = "$dir/31a.gtb";
#collect_models($f11, $f31, {1.1=>$d01, 1.2=>$d02, 1.5=>$d05, 1.9=>$d09});
my $f31 = "$dir/31.gtb";
#upgrade_crp_gtb($f31a, $f31);
my $f32 = "$dir/32_seq.gtb";
#gtb2Seq(-in=>$f31, -out=>$f32, -seq=>$f_seq, -opt=>1);


sub upgrade_crp_gtb {
    my ($fi, $fo) = @_;
    my $ti = readTable(-in=>$fi, -header=>1);
    for my $i (0..$ti->nofRow-1) {
        my ($srd, $locES, $locIS, $locCS, $loc5S, $loc3S) = map {$ti->elm($i, $_)} qw/strand locE locI locC loc5 loc3/;
        my ($locE) = parse_old_loc_str($locES);
        my ($locI) = parse_old_loc_str($locIS);
        my ($locC) = parse_old_loc_str($locCS);
        my ($loc5) = parse_old_loc_str($loc5S);
        my ($loc3) = parse_old_loc_str($loc3S);
        $srd = is_revsrd($srd, "+") ? "-" : "+";
        $ti->setElm($i, "strand", $srd);
        $ti->setElm($i, "locE", locAry2Str($locE));
        $ti->setElm($i, "locI", locAry2Str($locI));
        $ti->setElm($i, "locC", locAry2Str($locC));
        $ti->setElm($i, "loc5", locAry2Str($loc5));
        $ti->setElm($i, "loc3", locAry2Str($loc3));
        $ti->setElm($i, "note", $ti->elm($i, "cat3"));
        $ti->setElm($i, "cat3", $ti->elm($i, "cat2"));
        $ti->setElm($i, "cat2", "mRNA");
    }
    open(FH, ">$fo");
    print FH $ti->tsv(1);
    close FH;
}

__END__
my $f08_01 = file($d08, "01.tbl");
my $f09_01 = file($d09, "01.tbl");
#run R function assembly1()

my $f08_02 = file($d08, "02.tbl");
my $f09_02 = file($d09, "02.tbl");
#adjustCuration(-in1=>$f08_02, -in2=>$f09_01, -out=>$f09_02);


my $f05_22 = file($d05, "22_red.tbl");
my $f_gene_retire = file($DIR_Misc2, "jcvi", "09_genes_retire.txt");
#correct4RedGenes(-in1=>$f05_12, -in2=>$f_gene_retire, -out=>$f05_22);

