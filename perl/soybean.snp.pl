#!/usr/bin/perl
use strict; use Init; use Common; use Localdb; use Readfile; use Writefile;
use Bio::Seq; use Bio::SeqIO; use Path::Class;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
sub reformatGmGenome {
    my $fi  = file($DIR_Genome, "gm_101", "Gm_ori.fa");
    my $fo = file($DIR_Genome, "gm_101", "Gm.fa");
    my $seqhi = Bio::SeqIO->new(-file => $fi, -format => 'largefasta');
    my $seqho = Bio::SeqIO->new(-file => ">$fo", -format => "fasta");
    my $cnt = 0;
    while(my $seq = $seqhi->next_seq()) {
        print join("\t", $seq->id, $seq->description)."\n";
        my $seqNew = Bio::Seq->new(-id=>$seq->id, -seq=>$seq->seq);
        $seqho->write_seq($seqNew);
        $cnt ++;
    }
    print $cnt." sequences in $fi\n";
}
sub getSNPnum {
    #my ($fIn) = rearrange(['in'], @_);
    my $fIn = file($DIR_In, "Gm_markers.txt");
    my $fInH = new IO::File $fIn, "r";
    my $cntSnp = 0;
    while( <$fInH> ) {
        chomp;
        next if !$_;
        my @eleAry = split("\t", $_);
        $cntSnp ++ if $eleAry[1] eq "SNP";
    }
    print "$cntSnp seqs\n";
}
sub extract_snp_seq {
    my ($fIns, $fOut) = rearrange(['in', 'out'], @_);
    my $seqOH = Bio::SeqIO->new(-file => ">".$fOut, -format => "fasta");
    $fIns = [$fIns] if ref($fIns) eq "undefined";
    my $cntSeq = 0;
    for my $fIn (@$fIns) {
        my $seqIH = Bio::SeqIO->new(-file => "<".$fIn,  -format => "fasta");
        while(my $seq = $seqIH->next_seq()) {
            $cntSeq ++;
            $seq->id =~ /(ss\d+)/i;
            my $id = $1;
            die("cannot recognize ".$seq->id."\n") unless $id;
            $seq->description =~ /pos=(\d+).*len=(\d+).*subid="([\w\.\_\-]+)"/i;
            my ($pos, $len, $subid) = ($1, $2, $3);
            my $allele = "";
            $allele = $1 if $seq->description =~ /alleles="([A-Z\/\-]+)"/i;
            die("cannot extract info from ".$seq->description."\n") unless $subid;
            die("unequal length: $len <>".length($seq->seq)."\n") unless $len == length($seq->seq);
            my $desc = sprintf("%d/%d %s %s", $pos, $len, $subid, $allele);
            my $seqN = Bio::Seq->new(-id=>$id, -description=>$desc, -seq=>$seq->seq);
            $seqOH->write_seq($seqN);
        }
    }
}
sub get_snp_loc {
    my ($fi1, $fi2, $fo) = rearrange(['in1', 'in2', 'out'], @_);
    my $fhi = new IO::File $fi2, "r";
    my $fho = new IO::File $fo, "w";
    my $info = {};
    my $seqH = Bio::SeqIO->new(-file=>$fi1, -format=>'fasta');
    while(my $seq = $seqH->next_seq()) {
        my ($id, $desc) = ($seq->id, $seq->description);
        my @tmp = split(" ", $desc);
        my @locAry = split("/", $tmp[0]);
        $info->{$id} = [@locAry[0..1], @tmp[1..2]];
    }
    print scalar(keys %$info)." snps\n";
    print $fho join("\t", qw/id allele chr pos pct_idty alias note/)."\n";
    my $cnt = 0;
    while( <$fhi> ) {
        chomp;
        my $lh = blatParse($_);
        my $identity = ($lh->{matches}+1) / $lh->{qSize};
        my $id = $lh->{qName};
        my ($strand, $tName, $tStartAry, $qStartAry, $blockSizeAry) =
            ($lh->{strand}, $lh->{tName}, $lh->{tStartAry}, $lh->{qStartAry}, $lh->{blockSizeAry});
        die("$id not in seqInfo\n") unless exists $info->{$id};
        my ($pos, $qLen, $idL, $note) = @{$info->{$id}}[0..3];
        $pos = $qLen - $pos + 1 if $strand eq "-";
        my $hit = 0;
        my $posG;
        for my $i (0..$lh->{blockCount}-1) {
            my ($qStart, $tStart) = ($qStartAry->[$i], $tStartAry->[$i]);
            my $qEnd = $qStart + $blockSizeAry->[$i] - 1;
            if($pos >= $qStart && $pos <= $qEnd) {
                $hit = 1;
                $posG = $pos - $qStart + $tStart;
                last;
            } elsif($pos < $qStart) {
                my $posG = $tStart-1; #assumption
                last;
            }
        }
        $posG ||= $lh->{tEnd};
        my $flag = $hit ? "" : "Approx.";
        print $fho join("\t", $id, $note, $lh->{tName}, $posG, sprintf("%.02f", $identity), $idL, $flag)."\n";
    }
}
my $dir = dir($DIR_Misc3, "gm_snp1");
my $fis = [file($dir, "Q3_ss.fas"), file($dir, "Q4_ss.fas")];
my $f01 = file($dir, "01_seq.fa");
#extract_snp_seq(-in=>$fis, -out=>$f01);

#run gmap to locate SNPs on gm_101

my $f04 = file($dir, "04_snp_loc.txt");
#get_snp_loc(-in1=>$f01, -in2=>$f03, -out=>$f04);

$dir = dir($DIR_Misc2, "nbs", "gm");
my $f00 = file($dir, "nbs_lrr_nodule_Gm.gb");
my $f01 = file($dir, "01_na.fa");
my $f02 = file($dir, "02_aa.fa");
#nbs_gb_to_fasta($f00, $f01, $f02);
sub nbs_gb_to_fasta {
    my ($fi, $fo1, $fo2) = @_;
    my $seqhi = Bio::SeqIO->new(-file=>$fi, -format=>"genbank");
    my $seqho1 = Bio::SeqIO->new(-file=>">$fo1", -format=>"fasta");
    my $seqho2 = Bio::SeqIO->new(-file=>">$fo2", -format=>"fasta");
    while(my $seq = $seqhi->next_seq()) {
        my $id = $seq->accession_number;
        $id .= "_$1" if $seq->description =~ /cultivar\s+([\w\-]+)\s+/i;
        my $seq2 = Bio::Seq->new(-id=>$id, -seq=>$seq->seq);
        $seqho1->write_seq($seq2);
        $seqho2->write_seq($seq2->translate());
    }
}




