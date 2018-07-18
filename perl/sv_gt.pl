#!/usr/bin/perl -w
use strict; 
use FindBin;
use lib $FindBin::Bin;
use InitPath;
use Common;
use Seq;
use Data::Dumper;

my $f_refseq = "/project/youngn/zhoup/Data/genome/mt_35/41_genome.fa";
my $dir = "$DIR_Repo/mt_35/40_sv/41_shared";
my $f11 = "$dir/11.tbl";
my $d51 = "$dir/51_alt_refs";
my $f61 = "$dir/61_reads.tbl";
my $d62 = "$dir/62_reads";
my $d65 = "$dir/65_blat";
my $d66 = "$dir/66_exonerate";

#pindel2Tigrasv();
#prepare_alt_refs();
#print "bamPickBpReads -b $DIR_Repo/mt_35/40_sv/01_pos_sorted -p $f11 -o $f61\n";
#prepare_reads_seq();
#realign_reads();
sub prepare_alt_refs {
    my $t = readTable(-in=>$f11, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($id, $chr, $beg, $end, $beg_r, $end_r, $type, $size_d, $ins, $size_i, $n_ind) = $t->row($i);
        my $loc = Bio::Location::Simple->new(-start=>$beg-1000, -end=>$end+1000, -strand=>1, -seq_id=>$chr);
        my $seq_ref = seqRet($loc, $f_refseq);
        die "not $size_d + 2000: ".length($seq_ref)."\n" if length($seq_ref) != $size_d + 2000;
        my $seq_alt = substr($seq_ref, 0, 1000).$ins.substr($seq_ref, length($seq_ref)-1000, 1000);
        
        my $f_ref = file($d51, "$id.fa");
        my $seq1 = Bio::Seq->new(-id=>'ref', -seq=>$seq_ref);
        my $seq2 = Bio::Seq->new(-id=>'alt', -seq=>$seq_alt);
        writeSeq([$seq1, $seq2], $f_ref);
        print "$i done\n" if $i % 10 == 0;
    }
}
sub prepare_reads_seq {
    open(FHI, "<$f61");
    my $line = readline(FHI);
    my $seqHs;
    while(<FHI>) {
        chomp;
        my ($ind, $id_p, $type, $beg, $cigar, $id, $first, $seq) = split "\t";
        my $seq_read = Bio::Seq->new(-id=>"$ind|$id|$first", -seq=>$seq);
        my $seqH;
        if(!exists $seqHs->{$id_p}) {
            my $f_read = file($d62, "$id_p.fa");
            $seqH = Bio::SeqIO->new(-file=>">$f_read", -format=>'fasta');
            $seqHs->{$id_p} = $seqH;
        }
        $seqH = $seqHs->{$id_p};
        $seqH->write_seq($seq_read);
    }
}
sub realign_reads {
    my $t = readTable(-in=>$f11, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($id, $chr, $beg, $end, $beg_r, $end_r, $type, $size_d, $ins, $size_i, $n_ind) = $t->row($i);
        my $f_ref = file($d51, "$id.fa");
        die "cannot find refseq file for $id\n" unless -s $f_ref;
        my $f_reads = file($d62, "$id.fa");
        die "cannot find readseq file for $id\n" unless -s $f_reads;
        
        my $f_blat = file($d65, "$id.tbl");
#    run_pipe_blat($f_reads, $f_ref, $f_blat);
        my $f_exon = file($d66, "$id.tbl");
        validate_blat($f_reads, $f_ref, $f_blat, $f_exon, $size_d, $size_i);
#    last if $i == 0;
        print "$i done\n" if $i % 10 == 0;
    }
}
sub pindel2Tigrasv {
    my $f14 = file($dir, "14_tigrasv.tbl");
    open(FHO, ">$f14");
    my $t = readTable(-in=>$f11, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($id, $chr, $beg, $end, $beg_r, $end_r, $type, $size_d, $ins, $size_i, $n_ind, $n_reads, $n_reads_u, $indStr) = $t->row($i);
        my @ps = split(" ", $indStr);
        my @inds;
        for (@ps) {
            my @pps = split(":", $_);
            push @inds, $pps[0];
        }
        my $ind = join(",", @inds);
        if($type eq "D") {
            print FHO join("\t", $chr, $beg_r, $beg, $end, $end_r, "DEL", $size_d, "BWA", "SLX", $ind, "RP", $id)."\n";
        }
    }
}
sub run_pipe_blat {
    my ($f_reads, $f_ref, $fo) = @_;
    my $f_blat = "$DIR_Tmp/tmp.psl";
    my $f_sam = "$DIR_Tmp/tmp.sam";
    runCmd("blat $f_ref $f_reads $f_blat", 1);
    runCmd("psl2sam.pl $f_blat > $f_sam", 1);
    open(FH, "<$f_sam");
    my $h;
    while(<FH>) {
        chomp;
        my @ps = split "\t";
        my ($ridStr, $strand, $hit, $pos, $cigar, $score) = @ps[0,1,2,3,5,11];
        $score = substr($score, 5);
        die "unknown flag: $strand\n" if $strand != 0 && $strand != 16;
        $strand = $strand == 16 ? -1 : 1;
        $cigar =~ s/H/S/g;
        unless(exists $h->{$ridStr}->{$hit} && $h->{$ridStr}->{$hit}->[-1] > $score) {
            $h->{$ridStr}->{$hit} = [$pos, $strand, $cigar, $score];
        }
    }
    open(FHO, ">$fo");
    print FHO join("\t", qw/rid posR strandR cigarR scoreR posA strandA cigarA scoreA/)."\n";
    for my $rid (sort keys %$h) {
        my $h2 = $h->{$rid};
        my ($posR, $strandR, $cigarR, $scoreR, $mR) = ('', '', '', 0, 0);
        my ($posA, $strandA, $cigarA, $scoreA, $mA) = ('', '', '', 0, 0); 
        if(exists $h2->{ref}) {
            ($posR, $strandR, $cigarR, $scoreR) = @{$h2->{ref}};
        }
        if(exists $h2->{alt}) {
            ($posA, $strandA, $cigarA, $scoreA) = @{$h2->{alt}};
        }
        print FHO join("\t", $rid, $posR, $strandR, $cigarR, $scoreR, $posA, $strandA, $cigarA, $scoreA)."\n";
    }
}
sub validate_blat {
    my ($f_reads, $f_ref, $f_blat, $f_exon, $size_d, $size_i) = @_;
    open(FHO, ">$f_exon");
    print FHO join("\t", qw/acc rid first posR strandR cigarR scoreR posA strandA cigarA scoreA support/)."\n";
    my $t = readTable(-in=>$f_blat, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($ridStr, $posR, $strandR, $cigarR, $scoreR, $posA, $strandA, $cigarA, $scoreA) = $t->row($i);
        my ($cigarOpsR, $cigarOpsA) = map {cigarStr2Ops($_)} ($cigarR, $cigarA);
        my $mR = sum( map {($_->[1] =~ /[MI]/) ? $_->[0] : 0} @$cigarOpsR );
        my $mA = sum( map {($_->[1] =~ /[MI]/) ? $_->[0] : 0} @$cigarOpsA );
        my $dR = sum( map {($_->[1] =~ /[D]/) ? $_->[0] : 0} @$cigarOpsR );
        my $dA = sum( map {($_->[1] =~ /[D]/) ? $_->[0] : 0} @$cigarOpsA );
      
        my $tag = '';
        if($scoreR >= $scoreA && $scoreR >= 70) {
            if($posR > 1000-90+1 && $posR <= 1000+$size_d) {
                $tag = 'ref';
            } else {
                $tag = 'ref2';
            }
        } elsif($scoreR <  $scoreA && $scoreA >= 70) {
            if($posA > 1000-90+1 && $posA <= 1000+$size_i) {
                $tag = 'alt';
            }
        }
        my ($acc, $rid, $first) = split(/\|/, $ridStr);
        print FHO join("\t", $acc, $rid, $first, $posR, $strandR, $cigarR, $scoreR, $posA, $strandA, $cigarA, $scoreA, $tag)."\n";
    }
}
sub run_exonerate {
    my ($f_reads, $f_ref, $fo) = @_;
    my $cmd = "exonerate --model est2genome --showalignment F --showcigar T --showvulgar F --minintron 20 $f_reads $f_ref";
    my $h;
    open(OUT, "$cmd | ") or die "Failed: $!\n";
    while( <OUT> ) {
        chomp;
        next unless /^cigar/;
        my @ps = split " ";
        my ($qId, $qBeg, $qEnd, $qStrand, $hId, $hBeg, $hEnd, $hStrand, $score, @cigars) = @ps[1..$#ps];
        die "cigar not 2x\n" unless @cigars % 2 == 0;
        $hStrand = $hStrand eq "-" ? -1 : 1;
        if(exists($h->{$qId}->{$hId}) && $h->{$qId}->{$hId}->[-1] > $score) {
            #skip
        } else {
            my $pos = $hStrand eq "-" ? $hEnd+1 : $hBeg + 1;
            my @ops;
            push @ops, ["S", $qBeg-0] if $qBeg > 0;
            for my $i (0..(@cigars/2)-1) {
                push @ops, [$cigars[$i*2], $cigars[$i*2+1]];
            }
            push @ops, ["S", 90-$qEnd] if $qEnd < 90;
            @ops = reverse(@ops) if $hStrand == -1;
            $h->{$qId}->{$hId} = [$pos, $hStrand, \@ops, $score];
        }
    }
    open(FHO, ">$fo");
    print FHO join("\t", qw/rid posR strandR cigarR scoreR posA strandA cigarA scoreA/)."\n";
    for my $rid (sort keys %$h) {
        my $h2 = $h->{$rid};
        my ($posR, $strandR, $cigarR, $scoreR, $posA, $strandA, $cigarA, $scoreA) = ('', '', '', 0, '', '', '', 0);
        my ($opsR, $opsA) = ([], []);
        ($posR, $strandR, $opsR, $scoreR) = @{$h2->{ref}} if exists $h2->{ref};
        ($posA, $strandA, $opsA, $scoreA) = @{$h2->{alt}} if exists $h2->{alt};
        print FHO join("\t", $rid, $posR, $strandR, cigarOps2Str($opsR), $scoreR, $posA, $strandA, cigarOps2Str($opsA), $scoreA)."\n";
    }
}
sub cigarStr2Ops {
    my ($str) = @_;
    my @ops;
    while($str =~ /(\d+)([MIDNSHP])/ig) {
        push @ops, [$1, $2];
    }
    return \@ops;
}
sub cigarOps2Str {
    my ($ops) = @_;
    return @$ops ? join("", map {$_->[1].$_->[0]} @$ops) : "";
}



