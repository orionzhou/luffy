#!/usr/bin/perl -w
use strict;
use lib ($ENV{"m"}."/mt2/modules");
use Init;
use Common;
use Gff;
use Seq;
use Gtb;
use Genemodel;
use SignalP;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

$ENV{"SignalP"} = "/project/youngn/zhoup/Source/signalp-4.0";

my $f_seq = "$DIR_Misc3/spfinder/Athaliana/01_genome/01_refseq.fa";
my $f_gtb = "$DIR_Misc3/spfinder/Athaliana/01_genome/61_gene.gtb";
my $dir = "$DIR_Misc2/mapping/62_crp_at";

=cut
$f_seq = "$DIR_Misc3/spfinder/Osativa/01_genome/01_refseq.fa";
$f_gtb = "$DIR_Misc3/spfinder/Osativa/01_genome/61_gene.gtb";
$dir = "$DIR_Misc2/mapping/64_crp_os";

$f_seq = "$DIR_Misc3/spfinder/Mtruncatula/01_genome/01_refseq.fa";
$dir = "/project/youngn/zhoup/Data/misc2/crp/05_models";
=cut

my $f12 = "$dir/12.gtb";
my $f21 = "$dir/21.tbl";
#compare_models($f12, $f_gtb, $f21);

my $f01 = "$dir/01_seq.fa";
my $f22 = "$dir/22.tbl";
my $f31 = "$dir/31.gtb";
#combine_models(-fi=>$f22, -fo=>$f31, -gtb1=>$f12, -gtb2=>$f_gtb, -cat=>$f01);
my $f32 = "$dir/32_seq.gtb";
#gtb2Seq(-in=>$f31, -out=>$f32, -seq=>$f_seq, -opt=>1);

my $f41 = "$dir/41_sigp_score.tbl";
#sigp_score_gtb($f32, $f41);
my $f42 = "$dir/42_pep_score.tbl";
#pep_score_gtb($f32, $f42);
my $f45 = "$dir/45_stat.tbl";
#merge_stats($f32, $f41, $f42, $f45);

my $f51 = "$dir/51_true.gtb";
#pick_models($f32, $f45, $f51);
my $f52 = "$dir/52_true.gff";
#gtb2Gff($f51, $f52);



sub combine_models {
    my ($fi, $fo, $f_gtb1, $f_gtb2, $f_cat) = rearrange(['fi', 'fo', 'gtb1', 'gtb2', 'cat'], @_);
    my $ti = readTable(-in=>$fi, -header=>1);
    my $hg1 = readGtb(-in=>$f_gtb1, -opt=>2);
    my $hg2 = readGtb(-in=>$f_gtb2, -opt=>2);

    my $hc;
    open(FHI, "<$f_cat");
    while(<FHI>) {
        if(/^\>([\w\.\|\_]+) (CRP\d+)/) {
            $hc->{$1} = $2;
        }
    }
    close FHI;
    
    my $cnt = 1;
    open(FHO, ">$fo");
    print FHO join("\t", qw/id parent chr beg end strand locE locI locC loc5 loc3 phase source conf cat1 cat2 cat3 note/)."\n";
    for my $i (0..$ti->nofRow-1) {
        my ($id, $idG, $tag, $lenO, $len1, $len2, $note, $note2) = $ti->row($i);
        
        my $row;
        if($tag != 1 && $note eq "k") {
            die "no gene model for $id\n" unless exists $hg1->{$id};
            $row = $hg1->{$id};
            if($tag == 1) {
                my $rowG = $hg2->{$idG};
                $row->[0] = $rowG->[0]; 
                $row->[1] = $rowG->[1];
            }
            $row->[17] = "KS: $note2";
        } else {
            next if $tag == 9;
            die "no gene model for $idG\n" unless exists $hg2->{$idG};
            $row = $hg2->{$idG};
            $row->[17] = $note2 if $tag != 1;
        }
        
        $id =~ s/\_\d\.\d$//;
        if(!exists($hc->{$id})) {
            print "no cat for $id\n";
        } else {
            $row->[16] = $hc->{$id};
        }
        print FHO join("\t", @$row)."\n";
    }
    close FHO;
}
sub merge_stats {
    my ($f_gtb, $fs, $fp, $fo) = @_;

    my $tg = readTable(-in=>$f_gtb, -header=>1);
    my $ts = readTable(-in=>$fs, -header=>1);
    my $tp = readTable(-in=>$fp, -header=>1);

    open(FH, ">$fo");
    print FH join("\t", qw/id parent fam tag_sp score_sp n_cds lenC lenI codonStart codonStop preStop seq/)."\n";
    for my $i (0..$tg->nofRow-1) {
        my ($id, $pa, $fam, $seq) = map {$tg->elm($i, $_)} qw/id parent cat3 seq/;
        my ($idS, $tag_sp, $score_sp) = map {$ts->elm($i, $_)} qw/id tag score/;
        my ($idP, $codonStart, $codonStop, $preStop, $gap, $n_cds, $lenC, $lenI) = $tp->row($i);
        die "id conflict: $id != $idS [sigp]\n" unless $id eq $idS;
        die "id conflict: $id != $idP [pep]\n" unless $id eq $idP;
        print FH join("\t", $id, $pa, $fam, $tag_sp, $score_sp, $n_cds, $lenC, $lenI, $codonStart, $codonStop, $preStop, $seq)."\n";
    }
    close FH;
}
sub pick_models {
    my ($fi, $fs, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $ts = readTable(-in=>$fs, -header=>1);
    die "not equal no. of rows\n" if $t->nofRow != $ts->nofRow;

    my @idxs_rm;
    for my $i (0..$ts->nofRow-1) {
        my ($id, $pa, $fam, $tag_sp, $score_sp, $n_cds, $lenC, $lenI, $codonStart, $codonStop, $preStop, $seq) = $ts->row($i);
        my $cat1 = $t->elm($i, "cat1");
        push @idxs_rm, $i if !$tag_sp || $preStop || $cat1 eq "pseudogene"; 
    }

    printf "%d out of %d models removed\n", scalar(@idxs_rm), $t->nofRow;

    $t->delRows(\@idxs_rm);
    open(FH, ">$fo");
    print FH $t->csv(1, {delimiter=>"\t"});
    close FH;
}

