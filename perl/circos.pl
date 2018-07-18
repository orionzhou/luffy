#!/usr/bin/perl
use strict;
use FindBin;
use lib "$FindBin::Bin";

my ($feDb, $refDb) = ("mt_35", "mt_35");
getChrConf(-refdb=>$refDb, -opt=>1);
#getChrConf_2();
#linkPreProcess();
my $ps = {
      9 => {tag=>"29_all_loc", opt=>9, pres=>['25_nbs', '26_crp', '85_nod']},
    
    21 => {tag=>"31_loc"    , opt=>1, loc=> ["Mt5:8000000..9000000", "Mt2:31000000..34000000"] },
    22 => {tag=>"71_wgd_all", opt=>2, org=>"Mt", bylink=>1, format=>2},
    23 => {tag=>"72_wgd_ks" , opt=>2, org=>"Mt", bylink=>1, format=>2, ks=>[[0,0.2],[0.2,1.2],[1.2,3],[3,80]], cluster=>1_000_000},
    24 => {tag=>"81_crp_1"  , opt=>5, tag1=>'22', tag2=>"crp_1", tag=>"31"},
    25 => {tag=>"82_crp_2"  , opt=>5, tag1=>'22', tag2=>"crp_2", tag=>"32"},
    26 => {tag=>"83_seb"    , opt=>4, sp=>"mt", format=>2, cluster=>1_000_000},
    27 => {tag=>"84_wgd"    , opt=>3, ins=>2, sp=>"mt", format=>2, cluster=>1_000_000},
};
my $p = $ps->{9};
$p = { %$p, fedb=>$feDb, refdb=>$refDb };
#circosLoc($p);
#circosLink($p);
#circosBatch1($p);

#get_loc_nodup();
#get_loc_crp();
sub get_loc_nbs {
    my $ld = Localdb->new(-db=>$refDb);
    my $fi = file($DIR_Circos, "03_ids", "25_nbs.tbl");
    my $fo = file($DIR_Circos, "04_locs", "25_nbs.tbl");
    
    my $t = readTable(-in=>$fi, -header=>1);
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id chr beg end type type2/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $type) = $t->row($i);
        my $fe = $ld->getFeatureByName($id);
        print $fh join("\t", $id, $fe->seq_id, $fe->start, $fe->end, $type, "")."\n";
    }
}
sub get_loc_crp {
    my $fi = file($DIR_Misc2, "crp_plot/01_loc.txt");
    my $t = readTable(-in=>$fi, -header=>1);
    my $fo = file($DIR_Circos, "04_locs", "26_crp.tbl");
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id chr beg end type type2/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $c, $s, $e, $strand, $type2) = $t->row($i);
        my $type = $type2 ge "CRP1040" ? "NCR" : "nonNCR";
        print $fh join("\t", $id, $c, $s, $e, $type, $type2)."\n";
    }
}
sub get_loc_nodup {
    my $ld = Localdb->new(-db=>$refDb);
    my $fi = file($DIR_Circos, "03_ids", "85_nod.tbl");
    my $fo = file($DIR_Circos, "04_locs", "85_nod.tbl");
    
    my $t = readTable(-in=>$fi, -header=>1);
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id chr start end type/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id1, $id2, $nod1, $nod2) = $t->row($i);
        my ($cat1, $cat2);
        if( $nod1 eq "none" && $nod2 eq "none") {
            next;
        } elsif($nod1 ne "none" && $nod2 ne "none") {
            ($cat1, $cat2) = ("nodup_both") x 2;
        } elsif($nod1 eq "none" && $nod2 ne "none") {
            ($cat1, $cat2) = ("nodup_false", "nodup");
        } else {
            ($cat1, $cat2) = ("nodup", "nodup_false");
        }
        for my $i (0..1) {
            my $id = [$id1, $id2]->[$i];
            my $cat = [$cat1, $cat2]->[$i];
            my $fe = $ld->getFeatureByName($id);
            unless($fe) {
                $id =~ s/g/te/;
                $fe = $ld->getFeatureByName($id);
            }
            unless($fe) {print "$id cannot be found\n"; die;}
            print $fh join("\t", $id, $fe->seq_id, $fe->start, $fe->end, $cat)."\n";
        }
    }
}




