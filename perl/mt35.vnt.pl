#!/usr/bin/perl
use strict;
use FindBin;
use lib $FindBin::Bin;
use InitPath;
use Common; use Convert; use Localdb; 
use Readfile; use Writefile; use VntRet; use VntEffect;
use Bio::Seq; use Path::Class; use Data::Dumper;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
my $ps = {
  01 => { accs => [1..16, 19..21, 23..28, 101, 17..18, 22, 29, 30], 
    fedb => "mt_30", refdb => "mt_30", opt => {} },
  02 => { accs => [1..16, 19..21, 23..28, 101], 
    fedb => "mt_30", refdb => "mt_30", opt => {} },
  03 => { accs => [5, 6, 101, 29], 
    fedb => "mt_defl", refdb => "mt_30", opt => {} },
  04 => { accs => [5, 6, 101, 29], 
    fedb => "mt_30", refdb => "mt_30", opt => {} },
  05 => { accs => [1..16, 19..21, 23..28, 101], #deletion candidates
    fedb => "mt_35", refdb => "mt_35", opt => {} }, 
  06 => { accs => [1..16, 19..21, 23..28, 30..60, 101, 102, 323, 334], 
    fedb => "mt_35", refdb => "mt_35" },
  11 => { accs => [18, 29, 101], fedb => "mt_35", refdb => "mt_35", opt => {} },
};
my $opt = 06;
my $p = $ps->{$opt};
my ($feDb, $refDb) = map {$p->{$_}} qw/fedb refdb/;
$p->{accs} = [ map {sprintf "HM%03d", $_} @{$p->{accs}} ];

my $dirW = dir($DIR_Misc1, sprintf "seq%02d", $opt);
my $f01 = file($dirW, "01_id.tbl");
my $f02 = file($dirW, "02_loc.tbl");

my $vr = VntRet->new(-accs=>$p->{accs}, -covopt=>2, -fedb=>$feDb, -refdb=>$refDb);
#$vr->recSeqByIds(-in=>$f01, -out=>$dirW);

my $f31 = file($dirW, "31_vnt.tbl");
#$vr->getVntByIds(-in=>$f01, -out=>$f31);
my $f32 = file($dirW, "32_vnt_acc.tbl");
#splitVntPerAcc(-in=>$f31, -out=>$f32);
my $f35 = file($dirW, "35_vnt_sum.tbl");
#$vr->sumVnt(-in=>$f31, -out=>$f35, -opt=>2);

#$vr->getVntByLocs(-in=>$f02, -out=>$f31);
#$vr->sumVnt(-in=>$f31, -out=>$f35, -opt=>1);

my $ef = VntEffect->new(-fedb=>$p->{fedb}, -refdb=>$p->{refdb});
my $f33 = file($dirW, "36_vnt_effect.tbl");
#$ef->vntDesc(-in=>$f31, -out=>$f36);

my $f41 = file($dirW, "41_cov.tbl");
#$vr->getCovByIds(-in=>$f01, -out=>$f41);
#$vr->getCovByLocs(-in=>$f02, -out=>$f41);

my $f51 = file($dirW, '51_vntSum.txt');
#writeVntSum(-fid=>$f00, -fvnt=>$f21, -fcov=>$f41, -out=>$f51, -acc=>$accAry);



