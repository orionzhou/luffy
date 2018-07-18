#!/usr/bin/perl -w
BEGIN { push (@INC, $ENV{'mt'}."/mt2") };
use strict; use Init; use Common; use Localdb; use Run; 
use Annotate; use Readfile; use Writefile; use Seq;
use Path::Class; use Data::Dumper; use Parser; use VntWrite; use VntFilter;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
my $acc_31 = [ [1..16, 19..21, 23..28, 101], [17..18, 22, 29, 30] ];
my $acc_26 = [ [1..16, 19..21, 23..28, 101] ];
my $ps = {
  1 => { tag=>"run01", refdb=>'mt_30', win=>{size=>100_000, step=>100_000, snpsize=>1_500}, 
    co=>{maf=>0.1, gt=>0.75}, accsN=>$acc_26 },
  2 => { tag=>"run02", refdb=>'mt_30', win=>{size=>100_000, step=>100_000, snpsize=>1_500},
    co=>{maf=>0.1, gt=>0.75}, accsN=>$acc_26, cutoff=>{cov1=>2, cov2=>1, gt=>0.75} },
  3 => { tag=>"run03", refdb=>'mt_30', win=>{size=>100_000, step=> 50_000, snpsize=>1_500},
    co=>{maf=>0.1, gt=>0.75}, accsN=>$acc_26 },
  5 => { tag=>"run05", refdb=>'mt_30', win=>{size=>100_000, step=>100_000, snpsize=>1_500},
    co=>{maf=>0  , gt=>0   }, accsN=>$acc_26 },
};
my $p = $ps->{2};
$p->{accgroup} = makeAccAry($p->{accsN});
$p->{accs} = [mergeArray($p->{accgroup})];
my $dirW = dir($DIR_Misc1, $p->{tag});
my $f01 = file($dirW, "01_window.txt");
my $d02 = dir ($dirW, "02_simple_snp");
#getFilteredVnts(-do=>$d02, -fo=>$f01, -param=>$p);
#vntConvert(-dir=>$dirW, -format=>'fastphase');
sub getFilteredVnts {
  my ($fWin, $dirO, $p) = rearrange(['fo', 'do', 'param'], @_);
  system("mkdir -p $dirO") unless -d $dirO;
  my ($accG, $refDb, $co, $accs) = map {$p->{$_}} qw/accgroup refdb co accs/;
  my $fh = new IO::File $fWin, "w";
  print $fh join("\t", qw/round chr wStart wEnd cntSnp sStart sEnd/)."\n";
  for my $i (1..8) {
      my $chr = "chr$i";
      print "-> working on $chr\n";
      my $fi = dir($DIR_Variant, "07_snpcalled", "$chr.txt");
      my $vnt = readSsp($fi);
      $vnt = applyFilters(-vnt=>$vnt, -inds=>$accs, -co=>$co);
      my ($nInd, $nPos, $inds, $poss, $refs, $snps) = map {$vnt->{$_}} qw/nInd nPos inds poss refs snps/;
      my ($w, $idx) = assignWindows(-chr=>$chr, -param=>$p->{win}, -pos=>$poss, -refdb=>$refDb);
      for my $j (0..@$w-1) {
          my $w2 = $w->[$j];
          my ($wS, $wE, $swCnt) = @$w2;
          my $round = sprintf("%01d_%03d", $i, $j+1);
          if($swCnt == 0) {
              print $fh join("\t", $round, $chr, $wS, $wE, 0, "", "")."\n";
          }
          for my $k (0..$swCnt-1) {
              my $round2 = join("_", $round, $k+1);
              my @posIdxs = @{$idx->[$j]->[$k]};
              my $poss_win = [ @$poss[@posIdxs] ];
              my ($sS, $sE) = ($poss_win->[0], $poss_win->[-1]);
              my $vnt2 = filter_poss($vnt, $poss_win);
              my $pre = file($dirO, $round2);
              vntWrite(-pre=>$pre, -format=>'simplesnp', -vnt=>$vnt2);
              print $fh join("\t", $round2, $chr, $wS, $wE, scalar(@posIdxs), $sS, $sE)."\n";
          }
      }
  }
}
sub vntConvert {
  my ($dirW, $format) = rearrange(['dir', 'format'], @_);
  my $fi = file($dirW, "01_window.txt"); 
  my $t = readTable(-in=>$fi, -header=>1);
  my $h = { 'ldhat'=>11, 'ldhot'=>13, 'phase'=>21, 'fastphase'=>22, 'haploview'=>25, 'maxhap'=>31 };
  die "format not supported: $format\n" unless exists $h->{$format};
  my $dirO = dir($dirW, $h->{$format}."_$format", "01_in");
  system("mkdir -p $dirO") unless -s $dirO;
  for my $i (0..$t->nofRow-1) {
    my ($round, $chr, $wS, $wE, $cntSnp, $sS, $sE) = $t->row($i);
    next if $cntSnp < 2; 
    my $f1 = file($dirW, "02_simple_snp/$round.txt");
    my $ref = readSsp($f1);
    my $f2 = file($dirO, $round);
    vntWrite(-vnt=>$ref, -format=>$format, -pre=>$f2);
  }
}
my $d11 = dir($dirW, "11_ldhat");
#pipe_ldhat(-dir=>$d11, -fwin=>$f01);
my $d13 = dir($dirW, "13_ldhot");
#pipe_ldhot(-dir=>$d13, -fwin=>$f01);
my $d22 = dir($dirW, "22_fastphase" );
#pipe_fastphase(-dir=>$d22, -fwin=>$f01);
my $d25 = dir($dirW, "25_haploview" );
#pipe_haploview(-dir=>$d25, -fwin=>$f01);
my $d31 = dir($dirW, "31_maxhap" );
#pipe_maxhap(-dir=>$d31, -fwin=>$f01);




