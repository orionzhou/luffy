package Parser;
use strict;
use Common;
use Readfile;
use Path::Class; 
use Data::Table;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/parse_mummer_coords parse_fgenesh gsnap2sam gsnap2fa
    parse_fastphase parse_pindel
    filter_pindel/;
@EXPORT_OK = qw//;

sub parse_mummer_coords {
    my ($fi, $fo, $opt) = @_;
    die "$fi is not there\n" unless -s $fi;
    my $fhi = new IO::File $fi, "r";
    $opt ||= 1;
    my $cnt = 0;
    my @colNames = qw/s1 e1 s2 e2 strand1 strand2 len1 len2 pct_idty lenR lenQ covR covQ id1 id2/;
    my @cols;
    while( <$fhi> ) {
        chomp;
        next if ++$cnt <= 5;
        next unless $_;
        my @ps = split(/\|/);
        die @ps." cols [not 7] at\n\t".$_."\n" unless @ps == 7;
        map { $_ =~ s/(^\s+)|(\s+$)// } @ps;
        map { $_ =~ s/\s+/ / } @ps;
        my ($strand1, $strand2);
        my ($s1, $e1) = split(" ", $ps[0]);
        ($s1, $e1, $strand1) = $s1 <= $e1 ? ($s1, $e1, 1) : ($e1, $s1, -1);
        my ($s2, $e2) = split(" ", $ps[1]);
        ($s2, $e2, $strand2) = $s2 <= $e2 ? ($s2, $e2, 1) : ($e2, $s2, -1);
        my @tmp = ($s1, $e1, $s2, $e2, $strand1, $strand2);
        push @cols, [@tmp, map {split(" ", $_)} @ps[2..$#ps]];
    }
    my $t = new Data::Table(\@cols, \@colNames, 0);
    my $fho = new IO::File $fo, "w";
    if($opt == 1) {
        my $lens = $t->colRef("len1");
        my $pcts = $t->colRef("pct_idty");
        my @matches = map {int($lens->[$_] * $pcts->[$_] / 100)} (0..$t->nofRow-1);
        $t->addCol(\@matches, "matches", 0);
        my $cao = $t->colRefs([qw/matches pct_idty s1 e1 id1 s2 e2 id2/]);
        my $gan = [qw/score pct_idty s1 e1 id1 s2 e2 id2/];
        my $t2 = new Data::Table( $cao, $gan, 1 );
        $t = $t2;
        print $fho $t->tsv(0);
    } elsif($opt == 2) {
        my $newCol = [1..$t->nofRow];
        $t->addCol($newCol, "cnt", 0);
        $t->delCols([qw/lenR lenQ covR covQ/]);
        print $fho $t->tsv(0);
    }
}
sub parse_fgenesh {
    my ($fIn) = @_;
    my $fInH = new IO::File $fIn, "r";
    my $ref;
    while(<$fInH>) {
        chomp;
        if(/^\s*Seq name\:\s*(.+)$/) {
            ($ref->{id1}, $ref->{desc1}) = grep s/\s*$//, split(" ", $1, 2);
        } elsif(/^\s*Length of sequence\:\s*(\d+)/) {
            $ref->{len1} = $1;
        } elsif(/^\s*Homology\:\s*(.+)$/) {
            ($ref->{id2}, $ref->{desc2}) = grep s/\s*$//, split(" ", $1, 2);
        } elsif(/^\s*Length of homology\:\s*(\d+)/) {
            $ref->{len2} = $2;
        } elsif(/^\s*Number of predicted genes\s*(\d+)/) {
            $ref->{cntgene} = $1;
        } elsif(/^\s*no reliable predictions/) {
            $ref->{cntgene} = 0;
        } elsif(/^\s*(\d+)\s+([\+\-])\s+(\d+)\s+CDS/) {
            my @ps = split(" ");
            my ($idG, $strand, $type, $s, $e, $score, $s1, $e1) = @ps[0, 1, 3, 4, 6, 7, 8, 10];
            $ref->{gene}->{$idG} = {strand=>$strand, cds=>[]} unless exists $ref->{gene}->{$idG};
            my $refG = $ref->{gene}->{$idG};
            die "strand inconsistent\n".Dumper($ref) unless $strand eq $refG->{strand};
            my $refC = {start=>$s, end=>$e, score=>$score, orfstart=>$s1, orfend=>$e1, type=>$type};
            push @{$refG->{cds}}, $refC;
        }
    }
    return $ref;
}
sub get_sam_cigar {
    my ($ref) = @_;
    my @stats = @{$ref->{stat}};
    @stats = reverse @stats if $ref->{strand} == -1;
    my ($cigar, $mm, $tagindel) = ("", 0, 0);
    for my $stat (@stats) {
        my $match = $stat->{qe} - $stat->{qs} + 1;
        $cigar .= $match."M";
        $mm += $stat->{sub}; 
        if($stat->{del} > 0 || $stat->{ins} > 0) {
            if($tagindel == 0) {
                $cigar .= $stat->{ins}."I" if $stat->{ins} > 0;
                $cigar .= $stat->{del}."D" if $stat->{del} > 0;
                $tagindel = 1;
            } else {
                $tagindel = 0;
            }
        }
    }
    my $tag = "NM:i:$mm";
    return ($cigar, $tag);
}
sub get_sam_flag {
    my ($r1, $r2, $map1, $map2, $cnt) = @_;
    my ($flag1, $flag2) = (1 + 4*16, 1 + 8*16);
    if($map1 == 1) {
        die "map1=1, map2 != 1\n" unless $map2 == 1;
        $flag1 += 2;
    } elsif($map1 == 5) {
        $flag1 += 8;
    } elsif($map1 == 0) {
        $flag1 += 4;
    }
    if($map2 == 1) {
        $flag2 += 2;
    } elsif($map2 == 5) {
        $flag2 += 8;
    } elsif($map2 == 0) {
        $flag2 += 4;
    }
    if($r1 && $r1->{strand} == -1) {
        $flag1 += 1 * 16;
        $flag2 += 2 * 16;
    }
    if($r2 && $r2->{strand} == -1) {
        $flag1 += 2 * 16;
        $flag2 += 1 * 16;
    }
    if($cnt > 0) {
        $flag1 += 16*16;
        $flag2 += 16*16;
    }
    return ($flag1, $flag2);
}
sub make_sam_string {
    my ($id, $seq, $refs1, $refs2) = @_;
    my (@rst, $map1, $map2);
    my @refs1 = @$refs1;
    my @refs2 = @$refs2;
    if(@refs1) {
        my @tmp = uniq(map {$_->{tag}} @refs1);
        die "tags not unique: $id\n".Dumper(@refs1) unless @tmp == 1;
        $map1 = $tmp[0];
    } else {
        $map1 = 0;
    }
    if(@refs2) {
        my @tmp = uniq(map {$_->{tag}} @refs2);
        die "tags not unique: $id\n".Dumper(@refs2) unless @tmp == 1;
        $map2 = $tmp[0];
    } else {
        $map2 = 0;
    }
#  die "map1=$map1;map2=$map2\n".Dumper(@refs1).Dumper(@refs2) if $id eq "SNPSTER4:4:1:0:1170#0";
    if($map1 == 0) {
        my ($flag1, $flag2) = get_sam_flag(undef, undef, $map1, $map2, 0);
        my ($m, $ms) = @refs2 == 1 ? ($refs2[0]->{hid}, $refs2[0]->{stat}->[0]->{hs}) : ("*", 0);
        my $str1 = join("\t", $id, $flag1, "*", 0, 0, "*", $m, $ms, 0, $seq->{1}, "*");
        push @rst, $str1;
    } elsif($map1 == 5) {
            for my $i (0..$#refs1) {
                my $r1 = $refs1[$i];
                my ($flag1, $flag2) = get_sam_flag($r1, undef, $map1, $map2, $i);
                my ($cigar1, $tag1) = get_sam_cigar($r1);
                my $start1 = min(map {$_->{hs}} @{$r1->{stat}});
                my $seqStr1 = $r1->{strand} == -1 ? Bio::Seq->new(-seq=>$seq->{1})->revcom->seq : $seq->{1};
                my $str1 = join("\t", $id, $flag1, $r1->{hid}, $start1, 255, $cigar1, "*", 0, 0, $seqStr1, "*", $tag1);
                push @rst, $str1;
            }
    } else {
        die "not equal PE mappings: $id\n".Dumper(@refs1, @refs2) unless @refs1 == @refs2;
        for my $i (0..$#refs1) {
            my ($r1, $r2) = ($refs1[$i], $refs2[$i]);
            my ($flag1, $flag2) = get_sam_flag($r1, $r2, $map1, $map2, $i);
            my ($cigar1, $tag1) = get_sam_cigar($r1);
            my ($cigar2, $tag2) = get_sam_cigar($r2);
            my $start1 = min(map {$_->{hs}} @{$r1->{stat}});
            my $start2 = min(map {$_->{hs}} @{$r2->{stat}});
            my $seqStr1 = $r1->{strand} == -1 ? Bio::Seq->new(-seq=>$seq->{1})->revcom->seq : $seq->{1};
            my $seqStr2 = $r2->{strand} == -1 ? Bio::Seq->new(-seq=>$seq->{2})->revcom->seq : $seq->{2};
            my $str1 = join("\t", $id, $flag1, $r1->{hid}, $start1, 255, $cigar1, "=", $start2, $r1->{size}, $seqStr1, "*", $tag1);
            my $str2 = join("\t", $id, $flag2, $r2->{hid}, $start2, 255, $cigar2, "=", $start1, $r2->{size}, $seqStr2, "*", $tag2);
            push @rst, ($str1, $str2);
        }
    }
    if($map2 == 0) {
        my ($flag1, $flag2) = get_sam_flag(undef, undef, $map1, $map2, 0);
        my ($m, $ms) = @refs1 == 1 ? ($refs1[0]->{hid}, $refs1[0]->{stat}->[0]->{hs}) : ("*", 0);
        my $str2 = join("\t", $id, $flag2, "*", 0, 0, "*", $m, $ms, 0, $seq->{2}, "*");
        push @rst, $str2;
    } elsif($map2 == 5) {
        for my $i (0..$#refs2) {
            my $r2 = $refs2[$i];
            my ($flag1, $flag2) = get_sam_flag(undef, $r2, $map1, $map2, $i);
            my ($cigar2, $tag2) = get_sam_cigar($r2);
            my $start2 = min(map {$_->{hs}} @{$r2->{stat}});
            my $seqStr2 = $r2->{strand} == -1 ? Bio::Seq->new(-seq=>$seq->{2})->revcom->seq : $seq->{2};
            my $str2 = join("\t", $id, $flag2, $r2->{hid}, $start2, 255, $cigar2, "*", 0, 0, $seqStr2, "*", $tag2);
            push @rst, $str2;
        }
    }
    return join("\n", @rst)."\n";
}
sub parse_fastphase {
    my ($fi) = @_;
    my $fh = new IO::File $fi, "r";
    my @lines;
    my $tag = 0;
    while( <$fh> ) {
        chomp;
        $tag = 0 if /^END GENOTYPES/i;
        if($tag == 1 && $_) {
            push @lines, $_;
        }
        $tag = 1 if /^BEGIN GENOTYPES/i;
    }
    my $cntAcc = @lines / 3;
    die "unknown format $fi\n" unless $cntAcc =~ /^(\d+)$/;
    my $ref;
    for my $i (0..$cntAcc-1) {
        my ($acc, $ht1, $ht2) = @lines[$i*3..$i*3+2];
        my @a = split(" ", $ht1);
        my @b = split(" ", $ht2);
        die join("\n", $acc, $ht1, $ht2)."\n" unless @a == @b;
        $ref->{$acc} = [\@a, \@b];
    }
    return $ref;
}
sub parse_pindel {
    my ($fi, $fo) = @_;
    my $fhi = new IO::File $fi, "r";
    my $fho = new IO::File $fo, "w";
    print $fho join("\t", qw/id chr beg end beg_r end_r type size_d ins size_i n_ind n_reads n_reads_uniq ind/)."\n";
    my ($cnt, $cntI, $cntS) = (0, 0, 0);
    while(<$fhi>) {
        chomp;
        next unless /^\#/;
        my $line = readline($fhi);
        my @ps = split("\t", $line);
        my ($id, $typeStr, $ntStr, $chr, $b1, $e1, $b2, $e2, $n, $nu, $n_f, $n_fu, $n_r, $n_ru, $s1_str, $sum_str, $n_ind_all, $n_ind) = @ps[0..17];
        my ($type, $size_d) = split(" ", $typeStr);
        die "not a 'D(eletion)' type\n" if $type ne "D";
        my ($size_i, $ins) = ($1, $2) if $ntStr =~ /^NT (\d+) \"([ATCGN]*)\"/;
        $chr =~ s/ChrID //;
        $b1 =~ s/BP //;
        $b2 =~ s/BP\_range //;
        $b1 ++; $b2++; $e1 --; $e2 --;
        die "$id size inconsistent: $size_d\n" unless $size_d == $e1-$b1+1;
        $n =~ s/Supports //;
        $n_ind =~ s/NumSupSamples //;
        $n_f =~ s/\s//g;
        $n_r =~ s/\s//g;
        $cnt ++;
        
        my $h;
        for my $i (0..$n_ind_all-1) {
            my ($ind, $sn_f, $sn_fu, $sn_r, $sn_ru) = split(" ", $ps[19+$i]);
            my ($sn, $sn_u) = ($sn_f+$sn_r, $sn_fu+$sn_ru);
            next if $sn_u == 0;
            $h->{$ind} = [$sn, $sn_u];
        }
        my $ind_str = join(" ", map {sprintf "%s:%d:%d", $_, $h->{$_}->[0], $h->{$_}->[1]} sort(keys(%$h)));
        print $fho join("\t", $id, $chr, $b1, $e1, $b2, $e2, $type, $size_d, $ins, $size_i, $n_ind, $n, $nu, $ind_str)."\n";
        
        $line = readline($fhi);
        for my $j (1..$n) {
            $line = readline($fhi);
            my @ps = split(" ", $line);
            die "not 7 fileds at $id:\n$line\n" if @ps < 6;
            my ($strand, $matePos, $mapQ, $ind, $name1) = @ps[$#ps-4..$#ps];
            my ($name, $read) = split("/", substr($name1, 1));
#      my $beg = $b1 - length($str1);
#      my $cigar = sprintf("%dM%dD%dM", length($str1), $size_d, length($str2));
#      print $fho3 join("\t", $id, $ind, $name, $read, $beg, $cigar, $strand)."\n";
        }
    }
    printf "%d done\n", $cnt;
}



1;
__END__
