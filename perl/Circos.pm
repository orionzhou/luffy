package Circos;
use strict; use Init; use Common; use Localdb; use Run; use Annotate; use Readfile;
use Path::Class; use Data::Dumper; use Text::Template;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use CircosConf;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/; use Clone qw/clone/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/run_circos getChrConf getChrConf_2 
    chrStat recomRate
    getPlazaLegume linkPrProcess
    circosLoc circosLink
    clusterLinks 
    circosBatch1
    getOverlapLinks/;
@EXPORT_OK = qw//;
my @colorSet1 = qw/optblue optorange optyellow optpurple optgreen optred optviolet/;
my @colorSet2 = qw/nyt_blue nyt_red nyt_orange nyt_green nyt_yellow/;
sub run_circos {
    my ($fc) = @_;
    my $dirE = dir($DIR_Src, "circos-0.52");
    chdir $dirE || die "cannot cd to $dirE\n";
    my $cmd = "perl bin/circos --conf=$fc";
    runCmd($cmd, 1);
}
sub getChrConf {
    my ($refDb, $opt) = rearrange(['refdb', 'opt'], @_);
    $opt ||= 2;
    my $fo;
    if($opt == 1) {
        $fo = file($DIR_Circos, "01_conf", "assembly_$refDb.txt");
    } elsif($opt == 2) {
        $fo = file($DIR_R, "01_conf", "assembly_$refDb.txt");
    } else {
        die "unknown opt: $opt\n";
    }
    my $ld = Localdb->new(-db=>$refDb);
    my $fh = new IO::File $fo, "w";
    my $phaseColor = {0=>'ylgnbu-9-seq-5', 1=>"ylgnbu-9-seq-5", 2=>"ylgnbu-9-seq-5", 3=>"ylgnbu-9-seq-6"};
    print $fh join("\t", qw/id chr start end type/)."\n" if $opt == 2;
    my @fes = $ld->getFeatures(-types=>['chromosome']);
    for my $fe (sort {$a->id cmp $b->id} @fes) {
        my ($id, $s, $e) = ($fe->id, $fe->start - 1, $fe->end);
        print $fh join("\t", "chr", "-", $id, $id, $s, $e, $id)."\n" if $opt == 1;
    }
    @fes = $ld->getFeatures(-types=>['BAC', 'centromere']);
    for my $fe (sort {$a->seq_id cmp $b->seq_id || $a->start <=> $b->start} @fes) {
        my ($id, $chr, $start, $end, $type) = map {$fe->$_} qw/id seq_id start end primary_tag/;
        my ($color, $tag);
        if( $type eq 'centromere') {
            $color = 'black';
            $tag = $type;
        } elsif($type eq "BAC") {
            my $phase = $fe->score;
            $tag = "phase $phase";
            unless($phase) {
                $phase = 0;
                $tag = "?";
            }
            $color = $phaseColor->{$phase};
        }
        if($opt == 1) {
            print $fh join("\t", "band", $chr, $id, $id, $start-1, $end, $color)."\n";
        } elsif($opt == 2) {
            print $fh join("\t", $fe->id, $chr, $start, $end, $tag)."\n";
        }
    }
}
sub getChrConf_2 {
    my $di = dir($DIR_Circos, "02_coords");
    my $do = dir($DIR_Circos, "01_conf");
    my $fo = file($do, "rosids.txt");
    my $fh = new IO::File $fo, "w";
    my @sps = qw/mt gm lj vv/;
    for my $sp (@sps) {
        my $fi = file($di, "$sp.txt");
        my $t = readTable(-in=>$fi, -header=>1);
        my @chrs = uniq($t->col("chr"));
        $sp =~ s/^(\w)/\U$1/;
        for my $chr (@chrs) {
            next unless $chr =~ /\Q$sp\E\d+/i;
            my $t2 = $t->match_pattern("\$_->[5] eq '$chr'");
            my $end = max($t2->col("end"));
            print $fh join("\t", "chr", "-", $chr, $chr, 0, $end, $chr)."\n";
        }
    }
}
sub chrStat {
#my $statAry = [qw/snpDensity thetaPi thetaW tajimaD fuD/];
#for my $stat (@$statAry) {
#  my $fOut = file($DIR_Out, "stat_".$stat."_0.txt");
#  $fOutHs->{$stat} = new IO::File $fOut, "w";
#}
#chrStat(-in=>$fIn, -out=>$fOutHs);
    my ($fInHs, $fOutHs) = rearrange(['in', 'out'], @_);
    my @colNameAry = sort(keys(%$fOutHs));
    my $fOutTmp;
    for my $seq_id (sort(keys(%$fInHs))) {
        my $fInH = $fInHs->{$seq_id};
        while( <$fInH> ) {
            chomp;
            next if !$_ || $_ =~ /^region/;
            my @eleAry = split("\t", $_);
            my ($round, $start, $end, $length) = @eleAry[0..3];
            #print join("|", $seq_id, $start, $end, $length)."\n";
            my @stats = ($eleAry[4]/$length, $eleAry[5]/$length, $eleAry[6]/$length, $eleAry[8], $eleAry[9]);
            #print join("|", @stat)."\n";
            for my $i (0..@colNameAry-1) {
                die "no fileHandler to write\n" unless defined $fOutHs->{$colNameAry[$i]};
                $fOutTmp = $fOutHs->{$colNameAry[$i]};
                print $fOutTmp join("\t", $seq_id, $start, $end, $stats[$i])."\n";
            }
        }
    }
}
sub recomRate {
#my $fo = file($DIR_Out, 'recom.txt');
#recomRate(-in=>$fhis, -out=>$fo, -window=>100000);
    my ($fhis, $fo, $windowSize, $step) = rearrange(['in', 'out', 'window', 'step'], @_);
    my $fho = new IO::File $fo, "w";
    die "file not write-able: $fo" unless defined($fho);
    for my $seq_id (sort {substr($a, -1) <=> substr($b, -1)} keys %$fhis) {
        my $fhi = $fhis->{$seq_id};
        my $bacLoc = getWindows(-chr=>$seq_id, -window=>$windowSize, -step=>$step);
        my $pos = 0;
        for my $windowStart (sort {$a<=>$b} (keys(%$bacLoc))) {
            my $windowEnd = $bacLoc->{$windowStart};
            #print join(" - ", $seq_id, $windowStart, $windowEnd)."\n";
            my ($cntPos, $sumP, $sumD) = (0, 0, 0);
            my ($rateMean, $rateL, $rateU) = (0, 0, 0);
            while(my $line = readline($fhi)) {
                next if !$line || $line!~/^\d/;
                last if eof($fhi);
                my @eleAry = split("\t", $line);
                $pos = $eleAry[1];
                if($pos > $windowEnd) {
                    print $fho join("\t", $seq_id, $windowStart, $windowEnd, ($sumP*1000)/($sumD/1000))."\n" if $cntPos>0;
                    backOneLine($fhi);
                    last;
                } else {
                    #print join("\t", $pos, $windowEnd)."\n";
                    $cntPos ++;
                    $sumD += $eleAry[6];
                    $sumP += $eleAry[7];
                }
            }
        }
        #last if $seq_id eq "MtChr1";
    }
}
sub circosLoc {
    my ($p) = @_;
    my ($tag, $opt, $feDb, $refDb) = map {$p->{$_}} qw/tag opt fedb refdb/;
    my $fi = file($DIR_Circos, "03_ids",  "$tag.tbl");
    my $fo = file($DIR_Circos, "05_locs_circos",  "$tag.tbl");
    my $fh = new IO::File $fo, "w";
    my $ld = Localdb->new(-db=>$feDb) if $feDb;
    if($opt == 4) { #read from feDb::types features
        my $types = $p->{types};
        for my $type (@$types) {
            my @fes = $ld->getFeatures(-types=>$types);
            for my $fe (@fes) {
                print $fh join("\t", $fe->seq_id, $fe->start, $fe->end, $type)."\n";
            }
        }
    } elsif($opt == 9) { #merge locFile
        for my $pre (@{$p->{pres}}) {
            my $fi = file($DIR_Circos, "04_locs", "$pre.tbl");
            my $t = readTable(-in=>$fi, -header=>1);
            for my $i (0..$t->nofRow-1) {
                my ($id, $chr, $s, $e, $type1, $type2) = $t->row($i);
                my $type = $type2 ? "$type1:$type2" : $type1;
                $type =~ s/\s+/\_/g;
                print $fh join("\t", $chr, $s, $e, $type)."\n";
            }
        }
    } else {
        die "unknown opt: $opt\n";
    }
}
sub getMlpcInfo {
    my ($t) = @_;
    my $mid = $t->elm(0, "mid");
    my @cs1 = uniq($t->col("chr1"));
    my @cs2 = uniq($t->col("chr2"));
    die "mplc[$mid] has >1 chr1s\n" unless @cs1 == 1;
    die "mplc[$mid] has >1 chr2s\n" unless @cs2 == 1;
    my ($c1, $c2) = ($cs1[0], $cs2[0]);
#  $t->sort("start1", 0, 0);
#  my ($gs1, $ge1) = map {$t->elm($_, "gene1")} (0, $t->nofRow-1);
#  $t->sort("start2", 0, 0);
#  my ($gs2, $ge2) = map {$t->elm($_, "gene1")} (0, $t->nofRow-1);
#  print join("\t", $gs1, $ge1, $gs2, $ge2)."\n";
#  die "gene 1: $gs1 or $gs2\n" if $gs1 ne $gs2 && $gs1 ne $ge2;
#  die "gene 2: $ge1 or $ge2\n" if $ge1 ne $gs2 && $ge1 ne $ge2;
    my ($s1, $s2) = map { min($t->col($_)) } qw/start1 start2/;
    my ($e1, $e2) = map { max($t->col($_)) } qw/end1 end2/;
#  my $tr = $t->match_pattern("\$_->[10] <= 3");
#  my $ks = $tr->nofRow ? sum($tr->col("ks")) / $tr->nofRow : "";
#  $ks = sprintf("%.03f", $ks) if $ks;
    my $ks = join(" ", sort {$a<=>$b} $t->col("ks"));
    return ($mid, $c1, $s1, $e1, $c2, $s2, $e2, $ks);
}
sub linkPreProcess {
    my $dir = dir($DIR_Circos, "11_anchorpoints");
    my $fo1 = dir($dir, "01_links.txt");
    my $fo2 = dir($dir, "02_mlpcs.txt");
    my $fh1 = new IO::File $fo1, "w";
    my $fh2 = new IO::File $fo2, "w";
    my @headerL = qw/mid lid gene1 chr1 start1 end1 gene2 chr2 start2 end2 ks 4dtv/;
    my @headerM = qw/mid chr1 start1 end1 chr2 start2 end2 ks/;
    print $fh1 join("\t", @headerL)."\n";
    print $fh2 join("\t", @headerM)."\n";
    my @pres = qw/mt_mt gm_gm lj_lj vv_vv mt_gm mt_lj mt_vv gm_lj gm_vv lj_vv/;
    my $mid = 0;
    for( @pres ) {
        my $fi = file($dir, "$_.txt");
        next unless -s $fi;
        my $t1 = readTable(-in=>$fi, -header=>1);
        $t1->header([qw/mid genex geney ks 4dtv x y/]);
        my ($spx, $spy) = map {lc(substr($t1->elm(0, $_), 0, 2))} qw/genex geney/;
        my ($fx, $fy) = map {file($DIR_Circos, "02_coords", "$_.txt")} ($spx, $spy);
        my ($tx, $ty) = map {readTable(-in=>$_, -header=>1)} ($fx, $fy);
        my %hx = map {$tx->elm($_, "gene") => [$tx->elm($_, "chr"), $tx->elm($_, "start"), $tx->elm($_, "end")]} (0..$tx->nofRow-1); 
        my %hy = map {$ty->elm($_, "gene") => [$ty->elm($_, "chr"), $ty->elm($_, "start"), $ty->elm($_, "end")]} (0..$ty->nofRow-1); 
        my @midsL = uniq($t1->col("mid"));
        for my $midL (sort {$a<=>$b} @midsL) {
            $mid ++;
            my $t2 = $t1->match_pattern("\$_->[0] == $midL");
            my $tl = Data::Table->new([], \@headerL, 0);
            my @genexs = $t2->col("genex");
            my @geneys = $t2->col("geney");
            for my $i (0..$t2->nofRow-1) {
                my ($genex, $geney, $ks, $dtv) = map {$t2->elm($i, $_)} qw/genex geney ks 4dtv/;
                my $loc1 = $hx{uc $genex};
                my $loc2 = $hy{uc $geney};
                die "$genex not found\n" unless $loc1;
                die "$geney not found\n" unless $loc2;
                unless($ks =~ /^[\d\.]+$/) {
                    print "$genex $geney: $ks\n";
                    $ks = "";
                }
                $dtv = "" unless $dtv =~ /^[\d\.]+$/;
                $tl->addRow( [$mid, $i+1, $genex, @$loc1, $geney, @$loc2, $ks, $dtv] );
            }
            print $fh1 $tl->tsv(0);
            my @stats = getMlpcInfo($tl);
#      die join("\t", @stats)."\n" if $mid == 1;
            print $fh2 join("\t", @stats)."\n";
        }
        print "$fi done\n";
    }
    close $fh1;
    close $fh2;
}
sub extract_links {
    my ($p) = @_;
    my ($loc, $org, $ks, $dist, $refDb, $bylink) = map {$p->{$_}} qw/loc org ks dist refdb bylink/;
    my ($fi, $ic1, $is1, $ie1, $ic2, $is2, $ie2, $ik);
    if($bylink == 0) {
        $fi = file($DIR_Circos, "11_anchorpoints", "02_mlpcs.txt");
        ($ic1, $is1, $ie1) = (1..3);
        ($ic2, $is2, $ie2) = (4..6);
        $ik = 7;
    } elsif($bylink == 1) {
        $fi = file($DIR_Circos, "11_anchorpoints", "01_links.txt");
        ($ic1, $is1, $ie1) = (3..5);
        ($ic2, $is2, $ie2) = (7..9);
        $ik = 10;
    }
    my $t = readTable(-in=>$fi, -header=>1);
    if($org) {
        $t = $t->match_pattern("\$_->[$ic1] =~ /$org/i && \$_->[$ic2] =~ /$org/i");
        my $colO1 = $t->delCol("chr1");
        my @colN1 = grep s/$org/chr/i, @$colO1;
        $t->addCol(\@colN1, "chr1", $ic1);
        my $colO2 = $t->delCol("chr2");
        my @colN2 = grep s/$org/chr/i, @$colO2;
        $t->addCol(\@colN2, "chr2", $ic2);
    }
    if($loc) {
        my ($chr, $s, $e) = splitLocStr($loc);
        $t = $t->match_pattern("(\$_->[$ic1] eq '$chr' && \$_->[$is1] >= $s && \$_->[$ie1] <= $e) || (\$_->[$ic2] eq '$chr' && \$_->[$is2] >= $s && \$_->[$ie2] <= $e)");
    }
    if($ks) {
        my ($ksMin, $ksMax) = @$ks;
        $t = $t->match_pattern("\$_->[$ik] >= $ksMin && \$_->[$ik] < $ksMax");
    }
    if($dist) {
        my ($dMin, $dMax) = @$dist;
        $t = $t->match_pattern("\$_->[$ic1] eq \$_->[$ic2]");
        my $ld = Localdb->new(-db=>$refDb);
        my @idxD;
        my $h;
        for my $i (0..$t->nofRow-1) {
            my ($c, $s1, $s2) = map {$t->elm($i, $_)} qw/chr1 start1 start2/;
            unless(exists $h->{$c}) {
                my $loc = Bio::Location::Simple->new(-seq_id=>$c, -start=>1, -end=>getSeqLen($c, $refDb));
                my @fes = $ld->getFeatures(-types=>['gene', 'transposable_element_gene'], -loc=>$loc);
                my @starts = map {$_->start} @fes;
                @starts = sort {$a <=> $b} @starts;
                $h->{$c} = \@starts;
            }
            die "hash for $c not exist\n" unless exists $h->{$c};
            my $rank1 = bsearch($h->{$c}, $s1);
            my $rank2 = bsearch($h->{$c}, $s2);
            push @idxD, $i if abs($rank1-$rank2) > 100;
        }
        $t->delRows(\@idxD);
    }
    return $t;
}
sub extractLinksByLoc {
    my ($p) = @_;
    my ($loc, $tag) = map {$p->{$_}} qw/loc tag/;
    my $fi = file($DIR_Circos, "11_anchorpoints", "01_links.txt");
    my $t = readTable(-in=>$fi, -header=>1);
    my $fo = file($DIR_Circos, "12_links", $tag);
    my $tm = extract_links($p);
    my @mids = $tm->col("mid");
    my ($fo1, $fo2, $fo3) = map {"$fo\_$_.txt"} qw/links mlpcs conf/;
    my $fh1 = new IO::File $fo1, "w";
    my $fh2 = new IO::File $fo2, "w";
    my $fh3 = new IO::File $fo3, "w";
    my @mids_all = @mids;
    my $locs;
    for my $mid (@mids) {
        my $t2 = $t->match_pattern("\$_->[0] == $mid");
        my @stats;
        my ($s1, $s2) = map {floor($_/1_000_000)} @stats[1,4];
        my ($e1, $e2) = map {ceil ($_/1_000_000)} @stats[2,5];
        push @{$locs->{$stats[0]}}, [$s1, $e1];
        push @{$locs->{$stats[3]}}, [$s2, $e2];
        my $t4 = extract_link("$stats[0]:$stats[1]..$stats[2]", $t);
        push @mids_all, uniq($t4->col("mid"));
        my $t5 = extract_link("$stats[3]:$stats[4]..$stats[5]", $t);
        push @mids_all, uniq($t5->col("mid"));
    }
    for my $chr (keys %$locs) {
        for my $ref (@{posMerge($locs->{$chr})}) {
            printf $fh3 "<zoom>\nchr=%s\nstart=%du\nend=%du\nscale=5\n</zoom>\n",
                $chr, $ref->[0], $ref->[1];
        }
    } 
    for my $mid (uniq(@mids_all)) {
        my ($t3, @stats) = get_mlpc_info($mid, $t);
        print $fh2 join("\t", $mid, @stats[0..2])."\n";
        print $fh2 join("\t", $mid, @stats[3..5])."\n";
        for my $i (0..$t3->nofRow-1) {
            print $fh1 join("\t", "$mid\_$i", map {$t3->elm($i, $_)} qw/chr1 start1 end1/)."\n";
            print $fh1 join("\t", "$mid\_$i", map {$t3->elm($i, $_)} qw/chr2 start2 end2/)."\n";
        }
    }
    my @chrs = uniq($t->col("chr1"), $t->col("chr2"));
    print join(";", sort @chrs)."\n";
}
sub linkFromPlaza {
    my ($p) = @_;
    my ($pre, $bylink, $format) = map {$p->{$_}} qw/pre bylink format/;
    my $dirO = dir($DIR_Circos, "12_links");
    my $fo = file($dirO, "$pre.txt");
    my $fh = new IO::File $fo, "w";
    my $t = extract_links($p);
    if($bylink == 0) {
        $t = $t->subTable( [0..$t->nofRow-1], [qw/mid chr1 start1 end1 chr2 start2 end2 ks/] );
        $t->sort("chr1", 1, 0, "start1", 0, 0);
    } elsif($bylink == 1) {
        my $mids = $t->delCol( 'mid' );
        my $lids = $t->delCol( 'lid' );
        my @ids = map {join("_", $mids->[$_], $lids->[$_])} (0..$t->nofRow-1);
        $t->addCol(\@ids, "id", 0);
        $t = $t->subTable( [0..$t->nofRow-1], [qw/id chr1 start1 end1 chr2 start2 end2 ks/] );
    }
    print "\t".$t->nofRow." links extracted\n";
    $t = linkTrans($t, 2) if $format == 2;
    print $fh $t->tsv(0);
} 
sub linkFromGenePair {
    my ($p) = @_;
    my ($pre, $sp, $format, $opt) = map {$p->{$_}} qw/pre sp format opt/;
    my $fc = file($DIR_Circos, "02_coords", "$sp.txt");
    my $t = readTable(-in=>$fc, -header=>1);
    my %h = map {$t->elm($_, "gene") => [$t->elm($_, "chr"), $t->elm($_, "start"), $t->elm($_, "end")]} (0..$t->nofRow-1); 
    my $fi = file($DIR_Circos, "03_ids", "$pre.txt");
    my $fo = file($DIR_Circos, "12_links", "$pre.txt");
    my $t1 = readTable(-in=>$fi, -header=>0);
    my $to = Data::Table->new([], [qw/id chr start end value/], 0);
    for my $i (0..$t1->nofRow-1) {
        my ($id1, $id2, $value) = $t1->row($i);
        if($opt == 3) { 
            $id1 =~ s/medtr/mt/i;
            $id2 =~ s/medtr/mt/i;
        }
        next if $id1 =~ /^#/ || $id2 =~ /^#/;
        my $loc1 = $h{uc $id1};
        my $loc2 = $h{uc $id2};
        die "gene $id1 not found\n" unless $loc1;
        die "gene $id2 not found\n" unless $loc2;
        my ($c1, $s1, $e1) = @$loc1;
        my ($c2, $s2, $e2) = @$loc2;
        if($opt == 4) {
            die "$id1 and $id2 not on same chr\n" unless $c1 eq $c2;
        }
        $c1 =~ s/$sp/chr/i;
        $c2 =~ s/$sp/chr/i;
        if($opt == 3) { #normal gene pair format
            $to->addRow( [$i, $c1, $s1, $e1, "url=$value"] );
            $to->addRow( [$i, $c2, $s2, $e2, "url=$value"] );
        } elsif ($opt == 4) {  #seb's segdup block format
            $to->addRow([int($i/2)+1, $c1, $s1, $e2, "url=$value"]);
        } else {
            die "not opt 3 or 4: $opt\n";
        }
    }
    if($format == 1) {
        $to = linkTrans($to, 1);
        $to->sort("chr1", 1, 0, "start1", 0, 0);
    }
    my $fh = new IO::File $fo, "w";
    print $fh $to->tsv(0);
}
sub linkTrans {
    my ($ti, $opt) = @_;
    my $ncol = scalar( $ti->header() );
    my $h1 = [qw/id chr start end/];
    my $h2 = [qw/id chr1 start1 end1 chr2 start2 end2/]; 
    my $to;
    if($opt == 1) {
        die "not 4/5 cols\n" if $ncol != 4 && $ncol != 5;
        push @$h1, "value" if $ncol == 5;
        push @$h2, "value" if $ncol == 5;
        $ti->header($h1);
        $to = Data::Table->new([], $h2, 0);
        for my $i (0..$ti->nofRow/2-1) {
            my ($id1, $c1, $s1, $e1) = $ti->row($i * 2);
            my ($id2, $c2, $s2, $e2) = $ti->row($i * 2 + 1);
            die "$id1 != $id2 at row $i*2\n" unless $id1 eq $id2;
            my $row = [$id1, $c1, $s1, $e1, $c2, $s2, $e2];
            if($ncol == 5) {
                my $value = $ti->elm($i*2, 4);
                $value =~ s/^.*\=//;
                push @$row, $value;
            }
            $to->addRow($row);
        }
    } elsif($opt == 2) {
        die "not 7/8 cols\n" if $ncol != 7 && $ncol != 8;
        push @$h1, "value" if $ncol == 8;
        push @$h2, "value" if $ncol == 8;
        $ti->header($h2);
        $to = Data::Table->new([], $h1, 0);
        for my $i (0..$ti->nofRow-1) {
            my ($id, $c1, $s1, $e1, $c2, $s2, $e2) = $ti->row($i);
            my $row1 = [$id, $c1, $s1, $e1];
            my $row2 = [$id, $c2, $s2, $e2];
            push @$row1, "url=".$ti->elm($i, 7) if $ncol == 8;
            push @$row2, "url=".$ti->elm($i, 7) if $ncol == 8;
            $to->addRow($row1);
            $to->addRow($row2);
        }
    }
    return $to;
}
sub getOverlapLinks {
    my ($p) = @_;
    my ($tag, $tag1, $tag2) = map {$p->{$_}} qw/tag tag1 tag2/;
    my $fi1 = file($DIR_Circos, "12_links", "$tag1.txt");
    my $fi2 = file($DIR_Circos, "05_locs", "$tag2.txt");
    my $t1 = readTable(-in=>$fi1, -header=>0);
    $t1 = linkTrans($t1, 1);
    my $t2 = readTable(-in=>$fi2, -header=>0);
    $t2->header([qw/chr start end cat/]);
    my @idxs;
    for my $i (0..$t1->nofRow-1) {
        my ($id, $c1, $s1, $e1, $c2, $s2, $e2) = $t1->row($i);
#    my $t3 = $t2->match_pattern("\$_->[0] eq '$c1' && \$_->[1]>=$s1 && \$_->[2]<=$e1");
#    my $t4 = $t2->match_pattern("\$_->[0] eq '$c2' && \$_->[1]>=$s2 && \$_->[2]<=$e2");
        my $t3 = $t2->match_pattern("\$_->[0] eq '$c1' && ((\$_->[1]>=$s1 && \$_->[1]<=$e1) || (\$_->[2]>=$s1 && \$_->[2]<=$e1))");
        my $t4 = $t2->match_pattern("\$_->[0] eq '$c2' && ((\$_->[1]>=$s2 && \$_->[1]<=$e2) || (\$_->[2]>=$s2 && \$_->[2]<=$e2))");
        push @idxs, $i if $t3->nofRow > 0 && $t4->nofRow > 0;
    }
    $t1 = $t1->subTable(\@idxs, [$t1->header]);
    my $fo = file($DIR_Circos, "12_links", "$tag.txt");
    my $fh = new IO::File $fo, "w";
    print $fh linkTrans($t1, 2)->tsv(0);
}
sub clusterLinks {
    my ($p) = @_;
    my ($pre, $format, $cluster) = map {$p->{$_}} qw/pre format cluster/;
    my $dirI = dir($DIR_Circos, "12_links");
    my $dirO = dir($DIR_Circos, "13_links_merged");
    my $fi = file($dirI, "$pre.txt");
    my $fo = file($dirO, "$pre.txt");
    my $tl = readTable(-in=>$fi, -header=>0);
    my $tm = linkTrans($tl, 1);
    printf "  %4d links before clustering\n", $tm->nofRow;
    my @chrs = map {join("_", $tm->elm($_, "chr1"), $tm->elm($_, "chr2"))} (0..$tm->nofRow-1); 
    @chrs = uniq(@chrs);
    my $to = Data::Table->new([], [qw/id chr1 start1 end1 chr2 start2 end2/], 0);
    for (@chrs) {
        my ($c1, $c2) = split("_", $_);
        my $t2 = $tm->match_pattern("\$_->[1] eq '$c1' && \$_->[4] eq '$c2'");
        $t2->sort("start1", 0, 0);
        my @tags;
        my ($ps1, $pe1, $ps2, $pe2);
        for my $i (0..$t2->nofRow-1) {
            my ($s1, $e1, $s2, $e2) = map {$t2->elm($i, $_)} qw/start1 end1 start2 end2/;
            if($i == 0) {
                push @tags, 0;
                ($ps1, $pe1, $ps2, $pe2) = ($s1, $e1, $s2, $e2);
                next;
            }
            my $d1 = min( abs($ps1-$s1), abs($pe1-$e1), abs($ps1-$e1), abs($pe1-$s1) );
            my $d2 = min( abs($ps2-$s2), abs($pe2-$e2), abs($ps2-$e2), abs($pe2-$s2) );
            my $tag = ($d1 <= $cluster && $d2 <= $cluster) ? 1 : 0;
            ($ps1, $pe1, $ps2, $pe2) = ($s1, $e1, $s2, $e2);
            push @tags, $tag;
        }
        die "not ".$t2->nofRow." lines: ".@tags."\n" unless @tags == $t2->nofRow;
        my @idxRanges = getIdxRange(@tags);
        for (@idxRanges) {
            my ($is, $ie) = @$_;
            my $s1 = min( map {$t2->elm($_, "start1")} ($is, $ie) );
            my $s2 = min( map {$t2->elm($_, "start2")} ($is, $ie) );
            my $e1 = max( map {$t2->elm($_, "end1")} ($is, $ie) );
            my $e2 = max( map {$t2->elm($_, "end2")} ($is, $ie) );
            $to->addRow( [$t2->elm($is, "id"), $c1, $s1, $e1, $c2, $s2, $e2] );
        }
    }
    printf "  %4d links after clustering\n", $to->nofRow;
    $to = linkTrans($to, 2);
    my $fh = new IO::File $fo, "w";
    print $fh $to->tsv(0);
}
sub circosBatch1 {
    my ($p) = @_;
    my ($tag) = map {$p->{$_}} qw/tag/;
    my @pres = expandPres($p);
#  $cc->addPlot(-type=>'histogram', -file=>'jkj.txt', -r0=>'0.8r', -r1=>'0.9r', -min=>0, -max=>1, -color_fg=>"red", -color_bg=>"vlgrey");
    my $dirO = dir($DIR_Circos, "99_figs");
    for my $pre (@pres) {
        if($p->{cluster}) {
            my $fi1 = dir($DIR_Circos, "12_links", "$pre.txt");
            my $fi2 = dir($DIR_Circos, "13_links_merged", "$pre.txt");
            die "$fi1 is not there\n" unless -s $fi1;
            die "$fi2 is not there\n" unless -s $fi2;
            my $shows = [ [1, 0], [0, 1], [1, 1] ];
            for (@$shows) {
                my ($show1, $show2) = @$_;
                my $fp = "$dirO/$pre\_$show2$show1.png";
                my $cc = CircosConf->new(-out=>$fp);
                $cc->setChr(-optsp=>2);
                $cc->setLinks(-radius=>'0.8r');
                $cc->addLink(-tag=>"$pre\_1", -file=>$fi1, -ribbon=>1, -rules=>[[1, 'optorange', 1]]) if $show1 == 1;
                $cc->addLink(-tag=>"$pre\_2", -file=>$fi2, -ribbon=>1, -rules=>[[1, 'optgreen',  1]]) if $show2 == 1;
                my $fc = "/tmp/tmp.conf";
                $cc->write($fc);
                run_circos($fc);
            }
        } else {
            my $fi = dir($DIR_Circos, "12_links", "$pre.txt");
            my $fp = "$dirO/$pre.png";
            my $cc = CircosConf->new(-out=>$fp);
            $cc->setChr(-optsp=>2);
            $cc->setLinks(-radius=>'0.8r');
            $cc->addLink(-tag=>$pre, -file=>$fi, -ribbon=>1, -rules=>[[1, 'optblue', 1]]);
            my $fc = "/tmp/tmp.conf";
            $cc->write($fc);
            run_circos($fc);
        }
    }
}
sub expandPres {
    my ($p) = @_;
    my ($tag) = map {$p->{$_}} qw/tag/;
    my @pres;
    if(exists $p->{ks}) {
        for my $ks (@{$p->{ks}}) {
            push @pres, sprintf("%s_%.02f_%.02f", $tag, @$ks);
        }
    } elsif(exists $p->{ins}) {
        for my $i (1..$p->{ins}) {
            push @pres, sprintf("%s_%02d", $tag, $i);
        }
    } else {
        push @pres, $tag;
    }
    return @pres;
}
sub circosLink {
    my ($p) = @_;
    my ($tag, $opt) = map {$p->{$_}} qw/tag opt/;
    if($opt == 1) {
        extractLinksByLoc($p);
    } elsif($opt == 2) {
        my @ps;
        if( !exists $p->{ks} ) {
            my $p = clone($p);
            $p->{pre} = $tag;
            push @ps, $p;
        } else {
            for my $ks (@{$p->{ks}}) {
                my $p = clone($p);
                $p->{ks} = $ks;
                $p->{pre} = sprintf("%s_%.02f_%.02f", $tag, @$ks);
                push @ps, $p;
            }
        }
        for my $p (@ps) {
            my $ks = $p->{ks};
            print "Ks = ".join(" - ", @$ks)."\n";
            linkFromPlaza($p);
            clusterLinks($p) if exists $p->{cluster};
        }
    } elsif($opt == 3 || $opt == 4) {
        my @ps;
        if(!exists $p->{ins} || $p->{ins} == 0) {
            my $p = clone($p);
            $p->{pre} = $tag;
            push @ps, $p;
        } else {
            for my $i (1..$p->{ins}) {
                my $p = clone($p);
                $p->{pre} = sprintf("%s_%02d", $tag, $i);
                push @ps, $p;
            }
        }
        for my $p (@ps) {
            linkFromGenePair($p);
            clusterLinks($p) if exists $p->{cluster};
        }
    } elsif($opt == 5) {
        getLocFromFile($p);
    } else {
        die "unknonw opt in generating links: $opt\n";
    }
} 


1;
__END__
