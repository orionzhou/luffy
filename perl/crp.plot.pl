#!/usr/bin/perl -w
use strict;
use Cwd qw/abs_path/;
use File::Basename qw/dirname/;
BEGIN { unshift @INC, dirname(abs_path($0)); }
use File::Path qw/make_path remove_tree/;
use InitPath;
use Common;
use BioFeature;
use PlotFeature;

my $org = "Athaliana";
$org = "Mtruncatula_3.5";
#$org = "Osativa";
my $dir = "$DIR_Misc4/spada.crp.$org/41_perf_eval";
my $d41 = "$dir/41_imgs";
my @softs = qw/GeneID Augustus_de_novo GlimmerHMM GeneMark GeneWise_SplicePredictor Augustus_evidence SPADA All/;

my @ps = map {[$_, "$DIR_Misc4/spada.crp.$org/41_perf_eval/$_/61_final.gtb"]} @softs;
unshift @ps, ["At_Unannoated_SP_DB", "$DIR_Misc2/crp.ssp/$org/08.gtb"] if $org eq "Athaliana";
unshift @ps, ["OrysPSSP", "$DIR_Misc2/crp.ssp/$org/08.gtb"] if $org eq "Osativa";
unshift @ps, ["Current_Annotation", "$DIR_Misc4/spada.crp.$org/01_preprocessing/61_gene.gtb"];
unshift @ps, ["Gold_Stardard", "$DIR_Misc4/spada.crp.$org/41_perf_eval/01_model.gtb"];
@ps = map {[$_->[0] => readTable(-in=>$_->[1], -header=>1)]} @ps;

my @locs = extract_locs("$dir/21_stat_SPADA_0.001.tbl", "$dir/SPADA/61_final.tbl");
plot_by_locs(\@locs, $d41, \@ps);

sub extract_locs {
    my ($fs, $ft) = @_;

    my @ids;
    my $ts = readTable(-in=>$fs, -header=>1);
    for my $i (1..$ts->nofRow) {
        my ($id, $tag, $gene, $lenTP, $lenFP, $lenFN) = $ts->row($i-1);
        push @ids, $id if ($tag >= 3 && $tag <= 9) || ($tag == 2 && $lenFP+$lenFN>30);
        push @ids, $id;
    }
    my $hi = { map {$_=>1} @ids };

    my @locs;
    my $tt = readTable(-in=>$ft, -header=>1);
    for my $i (1..$tt->nofRow) {
        my ($id, $fam, $chr, $beg, $end) = $tt->row($i-1);
        push @locs, [$id, $chr, $beg, $end, $fam] if exists $hi->{$id};
    }
    printf " %d locs extracted\n", scalar(@locs);
    return @locs;
}
sub plot_by_locs {
    my ($locs, $do, $ps) = @_;
    make_path($do) unless -d $do;
#    remove_tree($do, {keep_root=>1});
    
    for my $i (1..@$locs) {
        my ($id, $chr, $beg, $end, $fam) = @{$locs->[$i-1]};
        $beg -= 1000;
        $end += 1000;
       
        my $fo = "$d41/$id.png";
        plot_transcripts_by_loc(-chr=>$chr, -beg=>$beg, -end=>$end, -out=>$fo, -p=>$ps);
        printf "  %5d / %5d\r", $i, scalar(@locs);
    }
    print "\n";
}


