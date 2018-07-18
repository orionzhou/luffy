#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use File::Path qw/make_path remove_tree/;
use InitPath;
use Mapping;
use Common;
use Data::Dumper;

my $ary = "ATH1";
my $db = "Athaliana";
#my $ary = "AtMtDEFL";
#my $db = "Mtruncatula_3.5";

my $dir = "$DIR_misc2/affy/$ary";
my $f05 = "$dir/05.fa";
#get_probeseq("$dir/01.tbl", $f05);

my $dir_db = "$DIR_data/db/blat";
my $f_db = "$dir_db/$db.2bit";
my $d11 = "$dir/11_mapping";
pipe_blat(-qry=>$f05, -tgt=>$f_db, -dir=>$d11, -tilesize=>8, -pctidty=>0.9, -best=>1);

sub get_probeseq {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $h;
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    for my $i (0..$t->nofRow-1) {
        my ($ids, $x, $y, $pos, $seq) = $t->row($i);
        $h->{$ids} ||= 0;
        $h->{$ids} ++;
        my $id = sprintf "%s.%02d", $ids, $h->{$ids};
        print FHO join("\n", ">$id", $seq)."\n";
    }
    close FHO;
}

