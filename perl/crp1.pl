#!/usr/bin/perl -w
use strict;
use Init;
use Common;
use Data::Dumper; 

my $f_seq = "/project/youngn/zhoup/Data/genome/mt_35/41_genome.fa";
my $f_gtb = "/project/youngn/zhoup/Data/genome/mt_35/10_model_Mt3.5v5/62_phase_fixed.gtb";

my $dir = "$DIR_Misc2/crp/03_hits";

my $f01a = "$DIR_Misc2/hmmsearch/01_crp_hmmsearchP/07.bes"; 
my $f01b = "$DIR_Misc2/hmmsearch/02_crp_hmmsearchX/07.bes";
my $f11 = "$dir/11.tbl";
my $f12 = "$dir/12_picked.tbl";
##hit_process($dir, {$f01a=>'hmmSearchP', $f01b=>'hmmSearchX'});

$dir = "$DIR_Misc2/crp/04_hits_picked";
my $f20 = "$dir/20_final.tbl";
my $f21 = "$dir/21_hits.tbl";
upgrade_loc($f20, $f21);

sub upgrade_loc {
    my ($fi, $fo) = @_;
    my $ti = readTable(-in=>$fi, -header=>1);
    open(FH, ">$fo");
    print FH join("\t", qw/id family chr beg end strand loc e source/)."\n";
    for my $i (0..$ti->nofRow-1) {
        my ($id1, $id2, $fam, $chr, $locS, $source, $e, $incons) = $ti->row($i);
        my ($beg, $end, $srd);
        my ($loc, $srd) = parse_old_loc_str($locS);
        $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
        $beg = $loc->[0]->[0];
        $end = $loc->[-1]->[1];
        print FH join("\t", $id1, $fam, $chr, $beg, $end, $srd, locAry2Str($loc), $e, $source)."\n";
    }
    close FH;
}

