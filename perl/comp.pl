#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.pl - pipeline of pairwise comparison between 2 genomes

=head1 SYNOPSIS
  
  comp.pl [-help] [-qry query-genome] [-tgt target-genome]

  Options:
    -h (--help)   brief help message
    -q (--qry)    query genome (def: HM056)
    -t (--tgt)    target genome (def: HM101)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use Cwd qw/abs_path/;
use List::Util qw/min max sum/;

my ($qry, $tgt) = ('PH207', 'Zmays_v4');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "qry|q=s" => \$qry,
  "tgt|t=s" => \$tgt,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $data = $ENV{'data'};
my $dirq = "$data/genome/$qry";
my $dirt = "$data/genome/$tgt";
my $qry_fas = "$dirq/11_genome.fas";
my $tgt_fas = "$dirt/11_genome.fas";
my $qry_2bit = "$dirq/21.blat/db.2bit";
my $tgt_2bit = "$dirt/21.blat/db.2bit";
my $qry_size = "$dirq/15.sizes";
my $tgt_size = "$dirt/15.sizes";
my $qry_size_bed = "$dirq/15.bed";
my $tgt_size_bed = "$dirt/15.bed";
my $qry_gap = "$dirq/16.gap.bed";
my $tgt_gap = "$dirt/16.gap.bed";

my $dir = "$data/misc3/$qry\_$tgt";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

prepare_blat();
##### qsub itasca
#process_blat1();
##### qsub itasca
#process_blat2();
#process_vnt();

sub prepare_blat {
  -d "01_seq" || make_path("01_seq");
  runCmd("breakseq.bygap.pl -i $qry_fas -o 00.fas -g 1000");
  runCmd("seq.splitlarge.py 00.fas 00.even.fas");
  runCmd("qsub.blat.pl -i 00.even.fas -o 01_seq -n 50 -t $tgt -g $qry");
}
sub process_blat1 {
  runCmd("cat 01_seq/part.*.psl > 11.psl");
  runCmd("psl2gal.pl -i 11.psl -o 11.gal");
  runCmd("gal.coord.pl -i 11.gal -p qry -q $qry_size -o - | \\
    gal.fix.ovlp.pl -i - -o - | \\
    gal.rmgap.pl -i - -q $qry_gap -t $tgt_gap -o - | \\
    gal.calib.pl -i - -q $qry_fas -t $tgt_fas -o 12.fixed.gal");
  runCmd("gal2psl.pl -i 12.fixed.gal -o 12.fixed.psl");

  runCmd("axtChain -linearGap=medium -psl 12.fixed.psl \\
    $tgt_2bit $qry_2bit 21.chain");
  runCmd("chainPreNet 21.chain $tgt_size $qry_size 23.chain");
  runCmd("chain2gal.pl -i 23.chain -o - | \\
    gal.calib.pl -i - -q $qry_fas -t $tgt_fas -o 23.gal");
  runCmd("gal2gax.pl -i 23.gal -o 23.gax");
  runCmd("gax2bed.pl -i 23.gax -p qry -o - | sortBed -i stdin | \\
    mergeBed -i stdin > 23.bed");
  runCmd("subtractBed -a $qry_size_bed -b $qry_gap | \\
    subtractBed -a stdin -b 23.bed | \\
    awk '(\$3-\$2) >= 50' - > 24.nov.bed");
  runCmd("seqret.pl -d $qry_fas -b 24.nov.bed -o 24.nov.fas");
  runCmd("rm 11.gal 23.chain 23.gax");
  runCmd("qsub.blat.pl -i 24.nov.fas -o 24.nov -n 1 -t $tgt -g $qry");
}

sub process_blat2 { 
  runCmd("cat 24.nov/part.*.psl > 24.nov.psl");
  runCmd("psl2gal.pl -i 24.nov.psl -o - | \\
    gal.coord.pl -i - -p qry -q $qry_size -o - | \\
    gal.fix.ovlp.pl -i - -o 25.gal");
  runCmd("gal2psl.pl -i 25.gal -o 25.psl");
  runCmd("cat 12.fixed.psl 25.psl > 31.1.psl");
  runCmd("pslSwap 31.1.psl 41.1.psl");
  runCmd("rm 25.gal");

#  pipe_pre1();
  pipe_chain_net("31", $qry_fas, $tgt_fas, $qry_2bit, $tgt_2bit, $qry_size, $tgt_size);
  pipe_chain_net("41", $tgt_fas, $qry_fas, $tgt_2bit, $qry_2bit, $tgt_size, $qry_size);
#  pipe_pre2();
  pipe_gal("31", $qry_fas, $tgt_fas, $qry_2bit, $tgt_2bit, $qry_size, $tgt_size);
  pipe_gal("41", $tgt_fas, $qry_fas, $tgt_2bit, $qry_2bit, $tgt_size, $qry_size);
}
sub pipe_pre1 {
  my $map = "$data/genome/$qry/raw.fix.fas.map"; 
  runCmd("rename.chain.pl -i ../23_blat.bak/31.2.chain \\
    -m $map -p qry -o 31.2.chain");
  runCmd("rename.chain.pl -i ../23_blat.bak/41.2.chain \\
    -m $map -p tgt -o 41.2.chain");
}
sub pipe_pre2 {
  my $map = "$data/genome/$qry/raw.fix.fas.map"; 
  runCmd("rename.pl -i ../23_blat.bak/31.3.gal -o - -m $map -c 7 -s 20 | \\
    tmp.addcol.pl -o 31.3.gal");
  runCmd("rename.pl -i ../23_blat.bak/31.5.gal -o - -m $map -c 7 -s 20 | \\
    tmp.addcol.pl | gal.addlev.pl -i - -n 31.5.net -o 31.5.gal");
  runCmd("rename.pl -i ../23_blat.bak/31.8.gal -o - -m $map -c 7 -s 20 | \\
    tmp.addcol.pl | gal.addlev.pl -i - -n 31.8.swap.net -o 31.8.gal");
  
  runCmd("rename.pl -i ../23_blat.bak/41.3.gal -o - -m $map -c 2 -s 20 | \\
    tmp.addcol.pl -o 41.3.gal");
  runCmd("rename.pl -i ../23_blat.bak/41.5.gal -o - -m $map -c 2 -s 20 | \\
    tmp.addcol.pl | gal.addlev.pl -i - -n 41.5.net -o 41.5.gal");
  runCmd("rename.pl -i ../23_blat.bak/41.8.gal -o - -m $map -c 2 -s 20 | \\
    tmp.addcol.pl | gal.addlev.pl -i - -n 41.8.swap.net -o 41.8.gal");
}
sub pipe_chain_net {
  my ($pre, $qFas, $tFas, $q2bit, $t2bit, $qSize, $tSize) = @_;
  runCmd("axtChain -linearGap=medium -psl $pre.1.psl \\
    $t2bit $q2bit $pre.2.chain");
  runCmd("chainPreNet $pre.2.chain $tSize $qSize $pre.3.chain");
  runCmd("chainSwap $pre.3.chain $pre.3.q.chain");

  runCmd("chainNet $pre.3.chain $tSize $qSize $pre.5.net $pre.5.q.net");
  runCmd("netChainSubset $pre.5.net $pre.3.chain stdout | \\
    chainSort stdin $pre.5.chain");
  runCmd("netChainSubset $pre.5.q.net $pre.3.q.chain stdout | \\
    chainSort stdin $pre.5.q.chain");

  runCmd("chainNet $pre.5.q.chain $qSize $tSize /dev/null $pre.8.net");
  runCmd("netChainSubset $pre.8.net $pre.3.chain $pre.8.chain");
}
sub pipe_gal {
  my ($pre, $qFas, $tFas, $q2bit, $t2bit, $qSize, $tSize) = @_;
  runCmd("chain2gal.pl -i $pre.5.chain -o - | \\
    gal.calib.pl -i - -q $qFas -t $tFas -o - | \\
    gal.addlev.pl -i - -n $pre.5.net -o $pre.5.gal");
  gal_expand("$pre.5.gal", "$pre.5", $qFas, $tFas, $qSize, $tSize);

  runCmd("chain2gal.pl -i $pre.8.chain -o - | \\
    gal.calib.pl -i - -q $qFas -t $tFas -o - | \\
    gal.addlev.pl -i - -n $pre.8.net -o $pre.8.gal");
  runCmd("gal.filter.pl -i $pre.8.gal -m 100 -p 0.6 -o $pre.9.gal");
  runCmd("gal2chain.pl -i $pre.9.gal -o $pre.9.chain");
  gal_expand("$pre.9.gal", "$pre.9", $qFas, $tFas, $qSize, $tSize);
}
sub gal_expand {
  my ($fi, $dir, $qFas, $tFas, $qSize, $tSize) = @_;
  $fi = abs_path($fi);
  -d $dir || make_path($dir);
  chdir($dir) || die "cannot chdir to $dir\n";

  runCmd("cp -f $fi gal");
  runCmd("gal.idx.pl -i gal -s $tSize");

  runCmd("gal2gax.pl -i gal -o gax");
  runCmd("gax.idx.pl -i gax -s $tSize");

  runCmd("gal2snp.pl -i gal -o snp -q $qFas -t $tFas");
  runCmd("snp.idx.pl -i snp -s $tSize");

  runCmd("gal2idm.pl -i gal -o idm");
  chdir "..";
}
sub process_vnt { 
  chdir "31.9" || die "cannot chdir to 31.9\n";
  runCmd("snp2vcf.pl -i snp -o snp.vcf -s $qry");
  runCmd("comp.sv.pl -q $qry -t $tgt");
  runCmd("vcf-concat snp.vcf ../../31_sv/11.sv.vcf | vcf-sort > vnt.1.vcf");
  runCmd("bcftools norm -c w -f $tgt_fas -o vnt.vcf vnt.1.vcf");
##  runCmd("vcf.fix.indel.pl -i vnt.2.vcf -o vnt.vcf");
  runCmd("vcf2tbl.pl -i vnt.vcf -o vnt.tbl");
  chdir "..";
} 

__END__

