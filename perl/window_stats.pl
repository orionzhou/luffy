#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use InitPath;
use Common;
use WindowStat;
use Path::Class;
my $acc_31 = [ [1..16, 19..21, 23..28, 101], [17..18, 22, 29, 30] ];
my $acc_26 = [ [1..16, 19..21, 23..28, 101] ];
my $ps = {
    1 => { tag=>"01_run02", refdb=>'mt_30', opt_gene=>1,  
        accsN=>$acc_26, cutoff=>{cov1=>2, cov2=>1, gt=>0.75} },
    2 => { tag=>"02_run02_candidates", refdb=>'mt_30', opt_gene=>1, 
        accsN=>$acc_26, cutoff=>{cov1=>2, cov2=>1, gt=>0.75} },
    3 => { tag=>"11_mt_35", refdb=>'mt_35', opt_gene=>2, 
        win=>{size=>100_000, step=>50_000} },
};
my $p = $ps->{3};
my $dirW = dir($DIR_Misc1, "stats_window", $p->{tag});
my $f01 = file($dirW, "01_window.tbl");
#write_windows(-out=>$f01, -p=>$p);
my $f02 = file($dirW, "02_gc.tbl");
#window_stats_gc(-fi=>$f01, -fo=>$f02, -refdb=>$p->{refdb});
my $f03 = file($dirW, "03_gene.tbl");
#window_stats_gene(-fi=>$f01, -fo=>$f03, -refdb=>$p->{refdb}, -opt=>$p->{opt_gene});
my $f04 = file($dirW, "04_cov.tbl");
#window_stats_cov(-fi=>$f01, -fo=>$f04, -p=>$p);
my $f05 = file($dirW, "05_gene_te.tbl");
#window_stats_features(-fi=>$f01, -fo=>$f05, -refDb=>$p->{refdb}, -types=>['gene', 'transposable_element_gene']);

my $f11_o = file($DIR_Misc1, "circos/04_locs/25_nbs.tbl");
my $f11 = file($dirW, "11_nbs.tbl");
#window_stats_loc(-fw=>$f01, -fi=>$f11_o, -fo=>$f11, -refDb=>$p->{refdb});
my $f12_o = file($DIR_Misc1, "circos/04_locs/26_crp.tbl");
my $f12 = file($dirW, "12_crp.tbl");
#window_stats_loc(-fw=>$f01, -fi=>$f12_o, -fo=>$f12, -refDb=>$p->{refdb});
$f11 = file($dirW, "11_genefam.tbl");
$f12 = file($dirW, "12_genefam_stat.tbl");
my $f_seq = file($DIR_Genome, $p->{refdb}, "41_genome.fa");
window_stats_loc(-fw=>$f01, -fi=>$f11, -fo=>$f12, -f_seq=>$f_seq);

my $f21 = file($dirW, "21_deletion.tbl");
my $f22 = file($dirW, "22_deletion.tbl");
#window_stats_loc(-fw=>$f01, -fi=>$f21, -fo=>$f22, -refDb=>$p->{refdb});

my $f51 = file($dirW, "51_stats.tbl");
#merge_stats(-ins=>[$f01, $f05, $f11, $f12], -out=>$f51);






