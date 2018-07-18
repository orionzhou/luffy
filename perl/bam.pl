#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Data::Dumper;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

use Bio::DB::Sam;

my $f_bam = "/home/youngn/zhoup/Data/misc2/mt_rnaseq/mt/accepted_hits.bam";
my $f_ref = "/home/youngn/zhoup/Data/genome/Mtruncatula_3.5/41_genome.fa";
my $sam = Bio::DB::Sam->new(-bam=>$f_bam, -fasta=>$f_ref);

my @tars = $sam->seq_ids;
my @alns = $sam->get_features_by_location(-seq_id => 'chr5',
    -start  => 7335001,
    -end    => 7343000);
for my $aln (@alns) {
    my ($beg, $end, $srd, $cigar) = map {$aln->$_} qw/start end strand cigar_str/;
    my $paired = $aln->get_tag_values('PAIRED');
    my $query_start = $aln->query->start;     
    my $query_end   = $aln->query->end;

    my $ref_dna   = $aln->dna;        # reference sequence bases
    my $query_dna = $aln->query->dna; # query sequence bases

    my @scores    = $aln->qscore;     # per-base quality scores
    my $match_qual= $aln->qual;       # quality of the match
    print join("\t", $beg, $end, $srd, $cigar)."\n";
}
