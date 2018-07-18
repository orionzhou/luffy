#!/usr/bin/perl
use strict;
use FindBin;
use lib "$FindBin::Bin";

use InitPath;
use Common;
use Run;
use Mtb;
use Hmm;
my $ps = {
  '02_crp_blastP' => {db=>"mt_35_v5_pro", program=>"blastp", e=>1e-3, pctidty=>0.5, best=>0},
};
my $tag = "02_crp_blastP";
my $p = $ps->{$tag};
my $dir = "$DIR_Misc2/blast/$tag";
system("mkdir -p $dir") unless -d $dir;

my $f_qry = "$dir/01_qry.fa"; 
my $f_ref = "$DIR_Genome/mt_35/41_genome.fa";
my $f_gtb = "$DIR_Genome/mt_35/10_model_Mt3.5v5/62_frame_fixed.gtb";

my $f02 = "$dir/02.bls";
#run_blast(-qry=>$f_qry, -out=>$f02, -param=>$p);
my $f03 = "$dir/03.mtb";
#bls2Mtb($f02, $f03);
my $f04 = "$dir/04_filtered.mtb";
#mtbFilter(-in=>$f03, -out=>$f04, -param=>$p);

my $f05 = "$dir/05_coord_recovered.mtb";
#recoverCoordP($f04, $f05, $f_gtb, $f_ref);
my $f06 = "$dir/06.til";
#hitTiling($f05, $f06);
my $f07 = "$dir/07.bes";
#hitCollapse($f06, $f07);

