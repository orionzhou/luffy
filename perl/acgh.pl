#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Math::Round qw/round/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $dirW = dir($DIR_Misc2, "acgh");
my $f00 = file($dirW, "00_raw.tbl");
my $f01 = file($dirW, "01_seq.fa");
#getProbeSeq($f00, $f01);
my $f02 = file($DIR_Misc2, "mapping/41_mt_35_cgh/06_all.mtb");
my $f05 = file($dirW, "05_mapping.tbl");
#filterProbeMappings($f02, $f05);
my $f06 = file($dirW, "06_mapping.bed");
#mtb2Bed(-in=>$f05, -out=>$f06);
print "cov_window -i $dirW/07_loc.tbl -o $dirW/08_cov.tbl -t opt23 -c 1\n";

print "grep '^chr5' \$data/genome/mt_35/10_model/63.bed > $dirW/51.bed\n";
#bed2Tbl(-in=>"$dirW/51.bed", -out=>"$dirW/51.tbl");
print "cov_window -i $dirW/51.tbl -o $dirW/52_cov.tbl -t opt23 -c 1\n";
print "intersectBed -wo -a $dirW/51.bed -b $dirW/06_mapping.bed > $dirW/54_intersect.bed\n";
my $f55 = file($dirW, "55_probe2gene.tbl");
#sum_bed_intersect("$dirW/54_intersect.bed", $f55);

my $acc = "HM029";
my $f10 = file($dirW, $acc, "10.tbl");
my $f21 = file($dirW, $acc, "21_sv_probes.tbl");
my $f22 = file($dirW, $acc, "22.tbl");
#pickIdentifiedProbes($acc, $f10, $f21, $f22);
my $f23 = file($dirW, $acc, "23_sv_gene.tbl");
#get_overlap_gene($f22, $f55, $f23);
my $f25 = file($dirW, $acc, "25_loc.tbl");
#pickProbeLocations($acc, $f05, $f22, $f25);
my $f31 = file($dirW, $acc, "31_cov.tbl");
#get_cov_stat($f15, $f31, $acc);

sub get_overlap_gene {
    my ($fi, $fg, $fo) = @_;
    my $fhg = new IO::File $fg, "r";
    my ($h1, $h2);
    while(<$fhg>) {
        chomp;
        my ($id, $chr, $beg, $end, $idG, $strand, $note, $len) = split("\t");
        die "$id mapped to >1 genes\n" if exists $h1->{$id};
        $h1->{$id} = [$idG, $len];
        $h2->{$idG} ||= [$chr, $beg, $end, $strand, $note];
    }
    
    my $t = readTable(-in=>$fi, -header=>1);
    my $fho = new IO::File $fo, "w";
    print $fho join("\t", qw/sv_id idG/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($sv_id, $n_probe, $probes) = map {$t->elm($i, $_)} qw/sv_id n_probe probes/;
        my @ids = split(" ", $probes);
        my $hg;
        for my $id (@ids) {
            next unless exists($h1->{$id});
            my ($idG, $len) = @{$h1->{$id}};
            $hg->{$idG} ||= [];
            push @{$hg->{$idG}}, [$id, $len];
        }
        my $hg2;
        for my $idG (keys %$hg) {
            my $rg = $h2->{$idG};
            my $lenG = $rg->[2] - $rg->[1];
            my $r = $hg->{$idG};
            my $len = sum( map {$_->[1]} @$r );
            my $str = join(" ", map {sprintf("%s[%d]", @$_)} @$r);
            if($len/$lenG > 0.3) {
                $hg2->{$idG} = [$len, $lenG];
            }
        }
        my $str = join(" ", keys %$hg2);
        print $fho join("\t", $sv_id, $str)."\n";
    }
}
sub bed2Tbl {
    my ($fi, $fo) = rearrange(['in', 'out'], @_);
    my $fhi = new IO::File $fi, "r";
    my $fho = new IO::File $fo, "w";
    print $fho join("\t", qw/chr beg end id strand/)."\n";
    while(<$fhi>) {
        chomp;
        my ($chr, $beg, $end, $id, $score, $strand, $note) = split("\t");
        $strand = $strand eq "-" ? -1 : 1;
        print $fho join("\t", $chr, $beg, $end, $id, $strand)."\n";
    }
}
sub sum_bed_intersect {
    my ($fi, $fo) = @_;
    my $fhi = new IO::File $fi, "r";
    my $h;
    while(<$fhi>) {
        chomp;
        next if /^\#/;
        my ($chr, $hBeg, $hEnd, $hId, $hStrand, $hNote, $chr2, $qBeg, $qEnd, $qId, $score, $qStrand, $len) = split("\t");
        my $stat = [$chr, $hBeg, $hEnd, $hId, $hStrand, $hNote, $len];
        if(exists($h->{$qId})) {
            my $stat_prev = $h->{$qId};
            if($stat_prev->[-1] > $stat->[-1]) {
                $stat = $stat_prev;
            }
        }
        $h->{$qId} = $stat;
    }
    my $fho = new IO::File $fo, "w";
    for my $id (keys %$h) {
        my @stats = @{$h->{$id}};
        print $fho join("\t", $id, @stats[0..6])."\n";
    }
}
sub mtb2Bed {
    my ($fi, $fo) = rearrange(['in', 'out'], @_);
    my $t = readTable(-in=>$fi, -header=>1);
    my $fh = new IO::File $fo, "w";
    print $fh "#track name=aCGH_mapping useScore=1\n";
    for my $i (0..$t->nofRow-1) {
        my ($qId, $qBeg, $qEnd, $strand, $hId, $hBeg, $hEnd, $matches, $qAlnLen, $qLen) = $t->row($i);
        my $pct_idty = sprintf "%.02f", $matches/$qAlnLen;
        $strand = $strand == -1 ? "-" : "+";
        print $fh join("\t", $hId, $hBeg-1, $hEnd, $qId, $pct_idty, $strand)."\n";
    }
}
sub getProbeSeq {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $seqh = Bio::SeqIO->new(-file=>">$fo", -format=>"fasta");
    for my $i (0..$t->nofRow-1) {
        my ($id, $seqStr) = map {$t->elm($i, $_)} qw/PROBE_ID PROBE_SEQUENCE/;
        $seqh->write_seq(Bio::Seq->new(-id=>$id, -seq=>$seqStr));
    }
}
sub filterProbeMappings {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $idx1 = first_index {$_ eq "hId"} $t->header;
    my $t2 = $t->match_pattern("\$_->[$idx1] =~ /^chr[1-8]\$/");
    
    my $ref = group($t2->colRef("qId"));
    my @rowsD;
    for my $qId (keys %$ref) {
        my ($idx, $nofId) = @{$ref->{$qId}};
        next if $nofId == 1;
        push @rowsD, ($idx..$idx+$nofId-1);
    }
    $t2->delRows(\@rowsD);
  
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/qId qBeg qEnd strand hId hBeg hEnd matches qAlnLen qLen/)."\n";
    my $cnt = 0;
    for my $i (0..$t2->nofRow-1) {
        my ($ma, $qId, $hId, $qLocStr, $hLocStr, $strand, $qGap, $hGap, $qAlnLen, $qLen) = 
            map {$t2->elm($i, $_)} qw/matches qId hId qLoc hLoc strand qGap hGap qAlnLen qLen/;
        my $qLoc = locStr2Obj($qLocStr, $qId);
        my $hLoc = locStr2Obj($hLocStr, $hId);
        my $tag = 1;
        $tag = 0 if $qLen - $ma > 0 || $qGap > 0 || $hGap > 0;
        if($tag == 1) {
            print $fh join("\t", $qId, $qLoc->start, $qLoc->end, $strand, $hId, $hLoc->start, $hLoc->end, $ma, $qAlnLen, $qLen)."\n";
            $cnt ++;
        }
    }
    close $fh;
    printf "%d are good\n", $cnt;
}
sub pickIdentifiedProbes {
    my ($acc, $f_ratio, $f_id, $f_out) = @_;
    my $tr = readTable(-in=>$f_ratio, -header=>1);
    my $ti = readTable(-in=>$f_id, -header=>1);
    my $fh = new IO::File $f_out, "w";
    print $fh join("\t", qw/sv_id n_probe cnv_len chr_o beg_o end_o probes ratios ratio_mean/)."\n";

    my $h_ratio = { map {$tr->elm($_, "PROBE_ID") => $tr->elm($_, "RATIO_CORRECTED")} (0..$tr->nofRow-1) };

    my $h = group($ti->colRef("Start"));
    my $cnv_cnt = 1;
    for my $key (sort {$a<=>$b} keys(%$h)) {
        my ($idx, $cnt) = @{$h->{$key}};
        my $cnv_id = sprintf "$acc\_%03d", $cnv_cnt++;
        my ($chr_o, $beg_o, $end_o) = map {$ti->elm($idx, $_)} qw/Chromosome Start End/;
        my $cnv_len = $end_o - $beg_o + 1;
        my @ids = map {$ti->elm($_, "Probe ID")} ($idx..$idx+$cnt-1);
        my @ratios = map {$h_ratio->{$_}} @ids;
        my $id_str = join(" ", @ids);
        my $ratio_str = join(" ", @ratios);
        my $mean_ratio = sum(@ratios) / @ratios;
        print $fh join("\t", $cnv_id, $#ids+1, $cnv_len, $chr_o, $beg_o, $end_o, $id_str, $ratio_str, $mean_ratio)."\n";
    }
}
sub pickProbeLocations {
    my ($acc, $f_mapping, $f_in, $f_out) = @_;
    my $tm = readTable(-in=>$f_mapping, -header=>1);
    my $ti = readTable(-in=>$f_in, -header=>1);
    my $fh = new IO::File $f_out, "w";
    print $fh join("\t", qw/id n_probe len len_o chr beg end/)."\n";

    my $h_loc = { map {$tm->elm($_, "qId") => [$tm->elm($_, "hBeg"), $tm->elm($_, "hEnd")]} (0..$tm->nofRow-1) };
    for my $i (0..$ti->nofRow-1) {
        my ($id, $n_probe, $len_o, $probes) = map {$ti->elm($i, $_)} qw/sv_id n_probe cnv_len probes/;
        my (@begs, @ends);
        for my $probe (split(" ", $probes)) {
            next unless exists $h_loc->{$probe};
            my $loc = $h_loc->{$probe};
            push @begs, $loc->[0];
            push @ends, $loc->[1];
        }
        my ($beg, $end) = (min(@begs), max(@ends));
        my $len = $end - $beg + 1;
        print $fh join("\t", $id, $n_probe, $len, $len_o, "chr5", $beg, $end)."\n";
    }
}
sub get_cov_stat {
    my ($f_loc, $f_cov, $acc) = @_;
    my $opt_ind;
    if($acc eq "HM018") {
        $opt_ind = "opt21";
    } elsif($acc eq "HM029") {
        $opt_ind = "opt22";
    } else {
        die "unsupported acc[$acc]\n";
    }
    runCmd("cov_window -i $f_loc -o $f_cov -t $opt_ind -c 1", 1);
}


