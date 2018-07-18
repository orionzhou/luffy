#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Path::Class; 
use Data::Dumper;

use InitPath;
use Common;
use Seq;

my $genome = "mt_35";
my $dir = dir($DIR_Genome, $genome);
my $f_seq = file($dir, "41_genome.fa");
my $f_mask = file($dir, "46_genome.mask.fa");
my $f_sum = file($dir, "51_sum.tbl");
#sum_genome($f_seq, $f_mask, $f_sum);
sub sum_genome {
    my ($f_seq, $f_mask, $f_sum) = @_;
    open(FH, ">$f_sum");
    print FH join("\t", qw/chr len_assembly len_bases length_alignable/)."\n";
    my $db = Bio::DB::Fasta->new($f_seq);
    my $db2 = Bio::DB::Fasta->new($f_mask);

    for my $id (sort $db->ids) {
        my $len_assembly = $db->length($id);
        my ($len_bases, $len_align) = (0) x 2;
        my $seq = $db->seq($id);
        while($seq =~ /([^N]+)/ig) {
            $len_bases += length($1);
        }

        my $len2 = $db2->length($id);
        die "length conflict: $len_assembly != $len2\n" unless $len_assembly == $len2;
        my $mask = $db2->seq($id);
        while($mask =~ /(0+)/g) {
            my $pos_end = pos($mask) - 1;
            my $pos_beg = $pos_end - length($1) + 1;
            $len_align += length($1);
        }
        print FH join("\t", $id, $len_assembly, $len_bases, $len_align)."\n";
    }
    close FH;
}

$dir = dir($DIR_Misc2, "genefam");
my $f_gtb = file($dir, "31_merged.gtb");
#sum_gene($f_gtb, "$dir/51_alignability.tbl");
sub sum_gene {
    my ($f_gtb, $fo) = @_;
    my $t = readTable(-in=>$f_gtb, -header=>1);
    open(FH, ">$fo");
    print FH join("\t", qw/id lenC lenM lenC_aln lenM_aln/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $parent, $chr, $beg, $end, $strand, $locC, $locE, $locU5, $locU3, $phase, $source, $conf, $cat1, $cat2, $cat3, $note) = $t->row($i);
        $locC = locStr2Obj($locC, $chr);
        $locE = locStr2Obj($locE, $chr);
        my $locM = Bio::Location::Simple->new(-seq_id=>$chr, -start=>$locE->start, -end=>$locE->end, -strand=>$locE->strand);
        my $seqC = maskRet($locC, $genome);
        my $seqM = maskRet($locM, $genome);

        my ($lenC, $lenM) = (length($seqC), length($seqM));
        my ($lenC_a, $lenM_a) = (0) x 2;
        while($seqC =~ /(0+)/g) { $lenC_a += length($1); }
        while($seqM =~ /(0+)/g) { $lenM_a += length($1); }
        print FH join("\t", $id, $lenC, $lenM, $lenC_a, $lenM_a)."\n";
    }
    close FH;
}


