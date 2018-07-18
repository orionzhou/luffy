#!/usr/bin/perl
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Data::Dumper; 
use InitPath;
use Common;
use Seq;
use Gff;
use Gtb;
use Bed;

sub processCpMt {
    my ($dir) = @_;
    my $f01 = "$dir/01_cp.gb";
    my $f02 = "$dir/02.fa";
    my $f03 = "$dir/03.gff";
#  gb2Gff(-in=>$f01, -out=>$f03, -outseq=>$f02, -type=>'chloroplase_sequence');
    my $f11 = "$dir/11_mt.gb";
    my $f12 = "$dir/12.fa";
    my $f13 = "$dir/13.gff";
#  gb2Gff(-in=>$f11, -out=>$f13, -outseq=>$f12, -type=>'mitochondrial_sequence');
    my $f41 = "$dir/41.fa";
    my $f51 = "$dir/51.gff";
}
sub merge_gene_genome {
    my ($f_seq, $f_gff, $fo) = @_;
    my $firstline = `head -n 1 $f_seq`;
    if($firstline eq "##FASTA\n") {
        system("cat $f_gff $f_seq > $fo");
    } else {
        system("echo '##FASTA' | cat $f_gff - $f_seq > $fo");
    }
    runCmd('sed -i \'s/\\t\w*_gene\\t/\\tgene\\t/g\' '.$fo);
}

my $org = "Mtruncatula_4.0";
my $dir = "$DIR_genome/$org";
my $f_seq = "$dir/11_genome.fa";

$dir = "$dir/gene";

my $f71 = "$dir/71.gff";

my $d80 = "$dir/80_cp_mt";
#processCpMt(-dir=>$d80);
my $f_seq_cpmt = "$d80/41.fa";
my $f_gff_cpmt = "$d80/51.gff";

my $f41 = "$dir/41_genome.fa";
#cat $f_seq_genome $f_seq_cpmt > $f41

my $f43 = "$dir/43_genome_snpEff.gff";
#merge_gene_genome($f41, $f_gff_gene, $f43);

