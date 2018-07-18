#!/usr/bin/perl
use strict;
use FindBin;
use lib "$FindBin::Bin";

use InitPath;
use Common;
use Bio::Seq;
use Data::Dumper;
use Seq;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use Time::HiRes qw/gettimeofday tv_interval/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
my $refDb = "mt_35";
my $dirW = dir($DIR_Misc3, "comparative");

my $chr = "chr3";
my $f01 = file($dirW, "01_$chr.fa");
#writeChrSeq($chr, $f01, $refDb);
#get_beds($chr);

my ($chr1, $chr2) = qw/chr5 chr3/;
my $f21 = file($dirW, "21_$chr1\_$chr2.align");
#mlpc2saf($chr1, $chr2, $f21);

sub writeChrSeq {
  my ($chr, $fo, $refDb) = @_;
  my $seqStr = seqRet(locStr2Obj("$chr:1..".getSeqLen($chr, $refDb)), $refDb);
#  $seqStr =~ s/[^ATCG]/N/ig;
  my $seq = Bio::Seq->new(-id=>$chr, -seq=>$seqStr);
  writeSeq(-seqs=>[$seq], -out=>$fo, -format=>'fasta');
}
sub mlpc2saf {
  my ($chr1, $chr2, $fo) = @_;
  my $f1 = file($DIR_Misc1, "circos/02_coords/mt.txt");
  my $t1 = readTable(-in=>$f1, -header=>1);
  my $h;
  for my $i (0..$t1->nofRow-1) {
    my ($idG, $idM, $beg, $end, $strand, $chr) = $t1->row($i);
    $strand = ($strand eq "-") ? -1 : 1;
    if($idG =~ /^MT(\d)G/) { $chr = "chr$1";};
    $h->{$idG} = [$chr, $beg, $end, $strand];
  }
  
  my $f2 = file($DIR_Misc1, "circos/11_anchorpoints/mt_mt.txt");
  my $t2 = readTable(-in=>$f2, -header=>1);
  my ($chrN1, $chrN2) = map {substr($_, length($_)-1, 1)} ($chr1, $chr2);
  $t2 = $t2->match_pattern("(\$_->[1] =~ /^MT$chrN1/ && \$_->[2] =~ /^MT$chrN2/) || (\$_->[1] =~ /^MT$chrN2/ && \$_->[2] =~ /^MT$chrN1/)");
  my $fh = new IO::File $fo, "w";
  for my $i (0..$t2->nofRow-1) {
    my ($mid, $id1, $id2, $ks) = $t2->row($i);
    my ($c1, $b1, $e1, $str1) = @{$h->{$id1}};
    my ($c2, $b2, $e2, $str2) = @{$h->{$id2}};
    if($c1 eq $chr2 && $c2 eq $chr1) {
      ($id1, $id2) = ($id2, $id1);
      ($c1, $c2) = ($c2, $c1);
      ($b1, $b2) = ($b2, $b1);
      ($e1, $e2) = ($e2, $e1);
      ($str2, $str1) = ($str1, $str2);
    }
    print $fh join("\t", $i+1, $b1, $e1, $b2, $e2, $str1, $str2, $id1, $id2, $ks, $mid)."\n";
  }
}
sub get_beds {
  my ($chr) = @_;
  runCmd("grep '^$chr' \$data/genome/$refDb/00_seq/32.bed > $dirW/10_bac_$chr.bed", 1);
  runCmd("grep '^$chr' \$data/genome/$refDb/10_model/63.bed > $dirW/11_gene_$chr.bed", 1);
  runCmd("grep '^$chr' \$data/misc2/igv/03_nbs.bed > $dirW/13_nbs_$chr.bed", 1);
  runCmd("grep '^$chr' \$data/misc2/igv/04_crp.bed > $dirW/14_crp_$chr.bed", 1);
}


