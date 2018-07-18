#!/usr/bin/perl -w
use strict; use Init; use Common; use Localdb; use Run; 
use Bio::Seq; use Bio::SeqIO; use Bio::SeqFeature::Generic;
use Readfile; use Writefile; use Align; use Parser; use Mapping;
use Gff; use Parser; use Seq; use Convert; use GeneModel;
use Time::HiRes qw/gettimeofday tv_interval/; use Data::Dumper; use Path::Class; 
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Math::Round qw/round/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
my $dir = dir($DIR_Misc3, "mt_mapping");
my ($refDb1, $refDb2) = qw/mt_35 mt_30/;
my $f01 = file($DIR_Misc2, "mapping/32_mt_35_map/01_seq.fa");
#writeGenomeSeq($refDb1, $f01);
my $f02 = file($DIR_Misc2, "mapping/32_mt_35_map/05_filtered.mtb");
my $f03 = file($dir, "03_mapping.tbl");
my $f04 = file($dir, "04_mapping.tbl");
#get_mapping($f02, $f03);
#recover_mapping_id($f01, $f03, $f04);

sub writeGenomeSeq {
    my ($refDb, $fo) = @_;
    my $seqh = Bio::SeqIO->new(-file=>">$fo", -format=>"fasta");
    my $cnt = 0;
    for my $chr (map {"chr".$_} (1..8)) {
        my $h = getWindows(-chr=>$chr, -winsize=>10000, -winstep=>10000, -db=>$refDb);
        for (@$h) {
            my ($beg, $end) = @$_;
            my $id = sprintf "%05d", ++$cnt;
            my $loc = Bio::Location::Simple->new(-seq_id=>$chr, -start=>$beg, -end=>$end, -strand=>1);
            my $locStr = "$chr:$beg-$end";
            my $seq = Bio::Seq->new(-id=>$id, -seq=>seqRet($loc, $refDb), -description=>$locStr);
            $seqh->write_seq($seq);
            print "$id\n";
        }
    }
}
sub get_mapping {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id qId qBeg qEnd strand hId hBeg hEnd/)."\n";

    my $cnt = 0;
    my $qIds = $t->colRef("qId");
    my $ref = group($qIds);
    for my $qId (keys %$ref) {
        my ($i, $nofId) = @{$ref->{$qId}};
        next if $nofId > 1;
        my ($id, $qId, $qLoc, $strand, $hId, $hLoc) = 
            map {$t->elm($i, $_)} qw/id qId qLoc strand hId hLoc/;
        my $qLocs = locObj2Ary(locStr2Obj($qLoc));
        my $hLocs = locObj2Ary(locStr2Obj($hLoc));
        $qLocs = [ sort {$a->[0] <=> $b->[0]} @$qLocs ]; 
        $hLocs = [ sort {$a->[0] <=> $b->[0]} @$hLocs ];
        $qLocs = [ reverse @$qLocs ] if $strand == -1;
        for my $i (0..@$qLocs-1) {
            my ($qb, $qe) = @{$qLocs->[$i]};
            my ($hb, $he) = @{$hLocs->[$i]};
            die "$id: $qb - $qe != $hb - $he\n" unless $qe - $qb == $he - $hb;

            print $fh join("\t", $id, $qId, $qb, $qe, $strand, $hId, $hb, $he)."\n";
            $cnt += $qe - $qb + 1;
        }
    }
    print "$cnt bps\n";
} 
sub recover_mapping_id {
    my ($fi1, $fi2, $fo) = @_;
    my $seqh = getSeqDesc($fi1);
    my $t = readTable(-in=>$fi2, -header=>1);
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/mt35_chr mt35_beg mt35_end strand mt30_chr mt30_beg mt30_end/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($aid, $qId, $qBeg, $qEnd, $strand, $chr2, $beg2, $end2) = $t->row($i);
        die "no seq desc for $qId\n" unless exists $seqh->{$qId};
        my $qDesc = $seqh->{$qId};
        my ($chr1, $beg1_g, $end1_g) = ($1, $2, $3) if $qDesc =~ /^(chr\d)\:(\d+)\-(\d+)$/;
        my $beg1 = $qBeg + $beg1_g - 1;
        my $end1 = $qEnd + $beg1_g - 1;
        print $fh join("\t", $chr1, $beg1, $end1, $strand, $chr2, $beg2, $end2)."\n";
    }
    close $fh;
}



