#!/usr/bin/perl
use strict;
use lib $ENV{"SCRIPT_HOME_PERL"};
use Data::Dumper;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use Time::HiRes qw/gettimeofday tv_interval/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

use InitPath;
use Common;
use CompareModel;

my $org;
$org = "Athaliana";
#$org = "Mtruncatula_3.5";
#$org = "Osativa";

my $dir = "$DIR_Misc2/crp.ssp/$org";
my $f_gtb_gs = "$DIR_Misc2/crp.gs/20130129/$org/05_renamed.gtb";

my $f_gtb_spa = "$DIR_Misc4/spada.crp.$org.simple/31_model_evaluation/61_final.gtb";
my $f_sta_spa = "$dir/11_stat_spa.tbl";
eval_gtb($f_gtb_gs, $f_gtb_spa, $f_sta_spa);

my $f_gtb_una = "$dir/08.gtb";
my $f_sta_una = "$dir/11_stat_una.tbl";
eval_gtb($f_gtb_gs, $f_gtb_una, $f_sta_una);

my $f_gtb_cur = "$DIR_Misc4/spada.crp.$org/01_preprocessing/61_gene.gtb";
my $f_sta_cur = "$dir/11_stat_cur.tbl";
eval_gtb($f_gtb_gs, $f_gtb_cur, $f_sta_cur);

sub model_eval {
    my ($f_eval, $f_ref_gtb, $fo) = @_;
    my $tv = readTable(-in=>$f_eval, -header=>1);
    my $tr = readTable(-in=>$f_ref_gtb, -header=>1);

    my $hg;
    for my $i (0..$tr->nofRow-1) {
        my ($id, $locS) = map {$tr->elm($i, $_)} qw/id locC/;
        my $loc = locStr2Ary($locS);
        my $len = locAryLen($loc);
        my $n_cds = @$loc;
        $hg->{$id} = [0, $len, $n_cds];
    }

    open(FH, ">$fo") || die "cannot open $fo for writing\n";
    print FH join("\t", qw/id tag gene lenTP lenFP lenFN exonTP exonFP exonFN/)."\n";
    for my $i (0..$tv->nofRow-1) {
        my ($id, $gene, $tag, $lenTP, $lenFP, $lenFN, $exonTP, $exonFP, $exonFN) = $tv->row($i);
        $tag = 5 if $tag == 7 || $tag == 8 || ($tag==2 && $lenFP+$lenFN>30);

        print FH join("\t", $id, $tag, $gene, $lenTP, $lenFP, $lenFN, $exonTP, $exonFP, $exonFN)."\n";
        $hg->{$gene}->[0] ++ if $gene ne "";
    }
    for my $gene (keys(%$hg)) {
        my ($cnt, $len, $n_cds) = @{$hg->{$gene}};
        if($cnt == 0) {
            print FH join("\t", '', 10, $gene, 0, 0, $len, 0, 0, $n_cds)."\n";
        } elsif($cnt > 1) {
            print "  $gene hit $cnt times\n";
        }
    }
    close FH;
}
sub eval_gtb {
    my ($f_gtb_qry, $f_gtb_ref, $f_stat) = @_;
    my $f_eval = "eval.tbl";
    compare_models($f_gtb_qry, $f_gtb_ref, $f_eval);
    model_eval($f_eval, $f_gtb_ref, $f_stat);
    runCmd("rm -f $f_eval");
}



