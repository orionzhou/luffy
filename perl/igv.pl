#!/usr/bin/perl -w
use strict;
use InitPath;
use Common; 
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $dir = "$DIR_Misc2/igv";
#tbl2Bed("$DIR_Misc2/nbs/mt_35/01_id.tbl", "$dir/03_nbs.bed");
#tbl2Bed("DIR_Misc2/crp_plot/01_loc.tbl", "$dir/04_crp.bed");
sub tbl2Bed {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    $t->sort("chr", 1, 0, "beg", 0, 0);
    
    open(FH, ">$fo") or die "cannot open $fo for reading\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $chr, $beg, $end, $strand, $fam) = $t->row($i);
        $strand = $strand == -1 ? "-" : "+";
        if($fam =~ /TIR/i) { $id = "$id\_$fam"; }
        print FH join("\t", $chr, $beg-1, $end, $id, 0, $strand)."\n"; 
    }
}

