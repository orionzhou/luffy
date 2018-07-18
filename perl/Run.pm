package Run;
use strict;
use Init;
use Common;
use Path::Class; 
use Data::Dumper;
use Seq;
use Time::HiRes qw/gettimeofday tv_interval/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/pipe_ldhat pipe_ldhot piep_fastphase pipe_haploview pipe_maxhap
    run_meme ldhot_varrec run_mummer/;
@EXPORT_OK = qw//;
sub pipe_ldhat {
    my ($dirW, $fwin) = rearrange(['dir', 'fwin'], @_);
    my $d01 = dir($dirW, "01_in");
    my $f03 = file($dirW, "03_convert.txt");
#  run_convert(-indir=>$d01, -fwin=>$fwin, -out=>$f03);
    my $f04 = file($dirW, "lk_26_t0.006.txt");
    my $d05 = dir( $dirW, "05_interval");
#  run_interval(-indir=>$d01, -outdir=>$d05, -flk=>$f04, -fwin=>$fwin);
    my $d06 = dir( $dirW, "06_stat");
#  run_stat(-indir=>$d05, -fwin=>$fwin, -outdir=>$d06);
    my $f07 = file($dirW, "07_recom");
    sum_stat(-indir1=>$d06, -indir2=>$d01, -out=>$f07, -fwin=>$fwin);
    my $d08 = dir( $dirW, "08_pairwise");
#  run_pairwise(-indir=>$d01, -outdir=>$d08, -flk=>$f04, -fwin=>$fwin);
}
sub run_interval {
    my ($dirI, $dirO, $fWin, $fi3) = 
        rearrange(["indir", "outdir", "fwin", "flk"], @_);
    die "$dirI is not there\n" unless -d $dirI;
    die "$fi3 is not there\n" unless -s $fi3;
    system("mkdir -p $dirO") unless -d $dirO;
    my $t = readTable(-in=>$fWin, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($round, $cntSnp) = map {$t->elm($i, $_)} qw/round cntSnp/;
        next unless $cntSnp > 1;
#    next unless $round =~ /^4/;
        my $fi1 = file($dirI, "$round.txt");
        my $fi2 = file($dirI, "$round\_loc.txt");
        die "$fi1 is not there\n" unless -s $fi1;
        die "$fi2 is not there\n" unless -s $fi2;
        my $pre = file($dirO, $round);
        print "$round\t$cntSnp SNPs\n";
        my $cmd = qq/interval -seq $fi1 -loc $fi2 -lk $fi3 -tag $pre -its 1000000 -bpen 5 -samp 2000/;
        runCmd($cmd);
        my $fo1 = file($dirO, "$round\_bounds.txt");
        my $fo2 = file($dirO, "$round\_rates.txt");
#    $round =~ s/\_1$//;
#    my $foo1 = file($dirO, "$round\_bounds.txt");
#    my $foo2 = file($dirO, "$round\_rates.txt");
#    die "$foo1 is not there\n" unless -s $foo1;
#    die "$foo2 is not there\n" unless -s $foo2;
#    system("mv $foo1 $fo1");
#    system("mv $foo2 $fo2");
    }
}
sub run_pairwise {
    my ($dirI, $dirO, $fWin, $fi3) = 
        rearrange(["indir", "outdir", "fwin", "flk"], @_);
    die "$dirI is not there\n" unless -d $dirI;
    die "$fi3 is not there\n" unless -s $fi3;
    system("mkdir -p $dirO") unless -d $dirO;
    my $t = readTable(-in=>$fWin, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($round, $cntSnp) = map {$t->elm($i, $_)} qw/round cntSnp/;
        next unless $cntSnp > 1;
        last if $i > 0;
        next unless $round =~ /^1/;
        my $fi1 = file($dirI, "$round.txt");
        my $fi2 = file($dirI, "$round\_loc.txt");
        die "$fi1 is not there\n" unless -s $fi1;
        die "$fi2 is not there\n" unless -s $fi2;
        my $pre = file($dirO, $round);
        my $t1 = [gettimeofday];
        print "$round\t$cntSnp SNPs\n";
        my $cmd = qq/pairwise -seq $fi1 -loc $fi2 -lk $fi3 -tag $pre/;
        open JJ, "| $cmd" || die "Failed: $! in \n$cmd\n";
        print JJ "1\n1000\n1001\n0\n1\n1\n0\n";
        while ( <JJ> ){
            chomp;
#      print $_."\n";
        }
        my $t2 = [gettimeofday];
        print "\t...done(".tv_interval($t1, $t2)."s)\n";
    }
}
sub run_stat {
    my ($dirI, $fWin, $dirO) = rearrange(["indir", 'fwin', "outdir"], @_);
    die "inDir is not there\n" unless -d $dirI;
    system("mkdir -p $dirO") unless -d $dirO;
    my $t = readTable(-in=>$fWin, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($round, $cntSnp) = map {$t->elm($i, $_)} qw/round cntSnp/;
        next unless $cntSnp > 1;
        my $fi  = file($dirI, "$round\_rates.txt");
        my $fo = file($dirO, "$round.txt");
        die "$fi is not there\n" unless -s $fi;
        my $pPath = file($DIR_Bin, "stat");
        my $cmd = qq/$pPath -input $fi -burn 150/;
        runCmd($cmd);
        move("res.txt", $fo);
        print "$round ...done\n";
    }
}
sub sum_stat {
    my ($dirI1, $dirI2, $fo, $fWin) = rearrange(['indir1', 'indir2', 'out', 'fwin'], @_); 
    die "$dirI1 is not there\n" unless -d $dirI1;
    die "$dirI2 is not there\n" unless -d $dirI2;
    my $fo1 = "$fo\_brief.tbl";
    my $fo2 = "$fo\_full.tbl";
    my $t = readTable(-in=>$fWin, -header=>1);
    my @keys = qw/1_mean 1_median 1_L95 1_U95 2_mean 2_median 2_L95 2_U95/;
    for (@keys) {
        $t->addCol([("") x $t->nofRow], $_);
    }
    my $fho2 = new IO::File $fo2, "w";
    print $fho2 join("\t", qw/round snp_id loc mean median L95 U95/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($round, $cntSnp, $s) = map {$t->elm($i, $_)} qw/round cntSnp sStart/;
        next unless $cntSnp > 1;
#    last if $round gt "1_004";
        my $fi1  = file($dirI1, "$round.txt");
        my ($stat, $t2) = sum_one_stat($fi1);
        for my $key (@keys) {
            my $value = $stat->{$key};
            $value /= ($t->elm($i, "sEnd") - $t->elm($i, "sStart") + 1) / 1_000 if $key =~ /^1/;
            $t->setElm($i, $key, sprintf("%.03f", $value));
        }
        my $fi2  = file($dirI2, "$round\_loc.txt");
        die "$fi2 is not there\n" unless -s $fi2;
        my $fhi2 = new IO::File $fi2, "r";
        my $tmp = readline($fhi2);
        my @loc;
        while( <$fhi2> ) {
            chomp;
            next unless $_;
            push @loc, $_ * 1000;
        }
        die "not ".@loc." loci: ".($t2->nofRow+1)." at $round\n" unless @loc == $t2->nofRow+1;
        for my $i (0..$#loc) {
            my @tmp = ("") x 4;
            @tmp = $t2->row($i-1) if $i > 0;
            print $fho2 join("\t", $round, $i+1, $s+$loc[$i]-1, @tmp)."\n";
        }
    }
    my $fho1 = new IO::File $fo1, "w";
    print $fho1 $t->tsv(1); 
}
sub sum_one_stat {
    my ($fi) = @_;
    die "$fi is not there\n" unless -s $fi;
    my $fhi = new IO::File $fi, "r";
    my $fTmp = "/tmp/stat.txt";
    my $fTH = new IO::File $fTmp, "w";
    my @keys = qw/loci mean median L95 U95/;
    while( <$fhi> ) {
        chomp;
        if(/^loci/i) {
            print $fTH join("\t", @keys)."\n";
        } else {
            print $fTH join("\t", split(" "))."\n";
        }
    }
    close $fTH;
    my $t = readTable(-in=>$fTmp, -header=>1);
    my $stat;
    for my $key (@keys) {
        $stat->{"1_$key"} = $t->elm(0, $key);
        my @tmp = $t->col($key);
        $stat->{"2_$key"} = sum(@tmp[1..$#tmp]) / $#tmp;
    }
    $t->delRow(0);
    $t->delCol("loci");
    system("rm $fTmp");
    return ($stat, $t);
}
sub run_convert {
    my ($dirI, $fWin, $fo) = rearrange(["indir", 'fwin', "out"], @_);
    die "$dirI is not there\n" unless -d $dirI;
    my $t = readTable(-in=>$fWin, -header=>1);
    my @statKeys = qw/segSites avgPwd thetaW varPwd tajimaD fuD/;
    for (@statKeys) {
        $t->addCol([("") x $t->nofRow], $_);
    }
    for my $i (0..$t->nofRow-1) {
        my ($round, $cntSnp) = map {$t->elm($i, $_)} qw/round cntSnp/;
        next unless $cntSnp > 1;
        my $fIn1 = file($dirI, "$round.txt");
        my $fIn2 = file($dirI, "$round\_loc.txt");
        my $stat = runConvert(-in1=>$fIn1, -in2=>$fIn2);
        print $round."\t";
        for my $statKey (@statKeys) {
            my $statValue = $stat->{$statKey};
            $t->setElm($i, $statKey, $statValue);
            die "cntSnp[$cntSnp] != segSites[$statValue]\n".Dumper($stat) if $statKey eq "segSites" && $cntSnp != $statValue;
        }
        print "done\n";
    }
    my $fh = new IO::File $fo, "w";
    print $fh $t->tsv(1);
}
sub runConvert { 
    my ($fIn1, $fIn2) = rearrange(["in1", "in2"], @_);
    die "$fIn1 is not there\n" unless -s $fIn1;
    die "$fIn2 is not there\n" unless -s $fIn2;
    open(JJ, "convert -seq $fIn1 -loc $fIn2 |") || die "Failed: $!\n";
    my $stat = {};
    while ( <JJ> ){
        chomp;
        if($_ =~ /Segregating sites.*?(\d+)/) {
            $stat->{segSites} = $1;
        } elsif($_ =~ /Average PWD.*?([\d\.]+)/) {
            $stat->{avgPwd}   = $1;
        } elsif($_ =~ /Watterson theta.*?([\d\.]+)/) {
            $stat->{thetaW}   = $1;
        } elsif($_ =~ /Tajima D statistic.*?([\d\.\-]+)/) {
            $stat->{tajimaD}  = $1;
        } elsif($_ =~ /Fu and Li D.*?([\d\.\-]+)/) {
            $stat->{fuD}      = $1;
        } elsif($_ =~ /Variance PWD.*?([\d\.]+)/) {
            $stat->{varPwd}   = $1;
        }
    }
    system("rm sites.txt");
    system("rm locs.txt");
    system("rm freqs.txt");
    return $stat;
}

sub pipe_ldhot {
    my ($dirW, $fwin) = rearrange(['dir', 'fwin'], @_);
    my $d01 = dir( $dirW, "01_in");
    my $d11_01 = dir($dirW, "../11_ldhat/01_in");
    my $f11_07 = dir($dirW, "../11_ldhat/07_recom_full.tbl");
    my $f02 = file($dirW, "02_conf.txt");
    my $d03 = dir( $dirW, '03_varrec');
#  prepare_ldhot(-indir1=>$d11_01, -indir2=>$d11_06, -out=>$d03, -fwin=>$fwin);
    prepare_ldhot(-in=>$f11_07, -out=>$d03, -fwin=>$fwin);
    my $f05 = file($dirW, '05_loc.txt');
    my $d06 =  dir($dirW, '06');
#  getRounds(-in1=>$f05, -in2=>$fwin);
#  run_ldhot(-indir1=>$d01, -indir2=>$d03, -outdir=>$d06, -fconfig=>$f02, -fwin=>$fwin, -rounds=>'3_126-4_001');
    my $f07 = file($dirW, "07_sum.txt");
#  sum_ldhot(-fwin=>$fwin, -indir=>$d06, -out=>$f07);
}
sub getRounds {
    my ($fi1, $fi2) = rearrange(['in1', 'in2'], @_);
    my $locH = readLocStr($fi1, "mt_30");
    my $t = readTable(-in=>$fi2, -header=>1);
    my $t2 = $t->match_pattern("\$_->[4] > 1");
    my @rst;
    for my $chr (keys %$locH) {
        my $ref1 = $locH->{$chr};
        for my $ref2 (@$ref1) {
            my ($start, $end) = @$ref2;
            my $t3 = $t2->match_pattern("\$_->[1] eq '$chr' && \$_->[5] >= $start && \$_->[6] <= $end");
            $t3->sort("round", 1, 0);
            printf "%s-%s\n", $t3->elm(0, "round"), $t3->elm($t3->nofRow-1, "round");
        }
    }
}
sub prepare_ldhot {
    my ($fi, $pre, $fWin) = rearrange(['in', 'out', 'fwin'], @_);
    my $dirO = $pre;
    my $fo1 = "$pre.tbl";
    my $fh1 = new IO::File $fo1, "w";
    print $fh1 join("\t", qw/round wStart wEnd cntSnp snp rho_mean rho_median rho_L95 rho_U95/)."\n";
    system("mkdir -p $dirO") unless -d $dirO;
    my $t = readTable(-in=>$fi, -header=>1);
    my $ref = group($t->colRef("round")); 
    for my $round (sort keys %$ref) {
        my ($idx_beg, $rows) = @{$ref->{$round}};
        next unless $rows > 1;
        my $idx_end = $idx_beg + $rows - 1;
        my $t2 = $t->subTable([$idx_beg..$idx_end], [$t->header]);
        my $ts = ldhot_window_sum($t2);

        my $fo2 = file($dirO, "$round.txt");
        my $fh2 = new IO::File $fo2, "w";
        for my $i (0..$ts->nofRow-1) {
            print $fh1 join("\t", $round, $ts->row($i))."\n"; 
            print $fh2 $ts->elm($i, "rho_mean")."\n";
        }
    }
}
sub ldhot_window_sum {
    my ($t) = @_;
    my $ts = Data::Table->new([], [qw/start end cntSnp snp rho_mean rho_median rho_L95 rho_U95/]);
    
    my @locs = $t->col("loc");
    my $winNumber = ceil( ($locs[-1] - $locs[0] + 1) / 1000 );
    my $h = [ map {[]} (1..$winNumber) ];
    for my $i (1..$#locs) {
        my $idx_win = ceil( ($locs[$i] - $locs[0] + 1) / 1000);
        push @{$h->[$idx_win-1]}, $i;
    }

    for my $i (0..@$h-1) {
        my $winS = $locs[0] + $i*1000;
        my $winE = $locs[0] + ($i+1)*1000 - 1;
        my @idxs = sort {$a<=>$b} @{$h->[$i]};

        my @stats;
        if(@idxs == 0) { #use next available snp
            my $idxWinNext = first_index {$_ > $i && @{$h->[$_]} > 0} (0..@$h-1);
            my $idx_next = $h->[$idxWinNext]->[0];
            @stats = map {$t->elm($idx_next, $_)} qw/mean median L95 U95/;
        } else {
            my ($posPrev, $len) = ($winS-1, 0);
            my ($sum_mean, $sum_median, $sum_l95, $sum_u95);
            for my $idx (@idxs) {
                my ($pos, $mean, $median, $l95, $u95) = map {$t->elm($idx, $_)} qw/loc mean median L95 U95/;
                my $len_inc = $pos - $posPrev;
                $len += $len_inc;
                $sum_mean   += $len_inc * $mean;
                $sum_median += $len_inc * $median;
                $sum_l95    += $len_inc * $l95;
                $sum_u95    += $len_inc * $u95;
                $posPrev = $pos;
            }
            if($idxs[-1] < $#locs) {
                my ($pos, $mean, $median, $l95, $u95) = map {$t->elm($idxs[-1]+1, $_)} qw/loc mean median L95 U95/;
                my $len_inc = $winE - $posPrev;
                $len += $len_inc;
                $sum_mean   += $len_inc * $mean;
                $sum_median += $len_inc * $median;
                $sum_l95    += $len_inc * $l95;
                $sum_u95    += $len_inc * $u95;
                die "$len != 1000\n" unless $len == 1000;
            }
            @stats = map {$_ / $len} ($sum_mean, $sum_median, $sum_l95, $sum_u95);
        }
        my @poss = map {$t->elm($_, "loc")} @idxs;
        @stats = map { sprintf "%.04f", $_ } @stats;
        $ts->addRow([$winS, $winE, $#idxs+1, join(" ", @poss), @stats]);
    }
    return $ts;
}
sub run_ldhot {
    my ($dirI1, $dirI2, $dirO, $fWin, $fi3, $rounds) = 
        rearrange(["indir1", "indir2", "outdir", "fwin", "fconfig", 'rounds'], @_);
    die "$dirI1 is not there\n" unless -d $dirI1;
    die "$dirI2 is not there\n" unless -d $dirI2;
    mkdir($dirO) unless -d $dirO;
    my $t = readTable(-in=>$fWin, -header=>1);
    my ($rs, $re) = split("-", $rounds);
    my $t2 = $t->match_pattern("\$_->[0] ge '$rs' && \$_->[0] le '$re'");
    for my $i (0..$t2->nofRow-1) {
        my ($round, $cntSnp) = map {$t2->elm($i, $_)} qw/round cntSnp/;
        next if $cntSnp < 2;
        my $fi1 = file($dirI1, "$round\_ldhot.txt");
        my $fi2 = file($dirI2, "$round\_varrec.txt");
        die "$fi1 is not there\n" unless -s $fi1;
        die "$fi2 is not there\n" unless -s $fi2;
        my $fo = file($dirO, "$round.txt");
        my $ti1 = [gettimeofday];
        print "$round\t$cntSnp SNPs\n";
        my $cmd = join(" ", "ldhot", $fi3, $fi1, "-V", $fi2, $fo);
        runCmd($cmd, 0);
        if( -s "$fo.sum" ) {
            system("mv $fo.sum $fo");
            my $ti2 = [gettimeofday];
            print "\n\t...done(".tv_interval($ti1, $ti2)."s)\n";
        } else {
            print "\tround[$round] failed\n";
        }
    }
}
sub sum_ldhot {
    my ($fw, $dirI, $fo) = rearrange(['fwin', 'indir', 'out'], @_);
    my $t = readTable(-in=>$fw, -header=>1);
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/chr start end LR rhohat/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($round, $cntSnp, $chr, $sS, $sE) = map {$t->elm($i, $_)} qw/round cntSnp chr sStart sEnd/; 
        next if $cntSnp < 2;
        next unless $round =~ /^[235]/;
        my $fi = file($dirI, "$round.txt");
        if( -s $fi ) {
            my $fhi = new IO::File $fi, "r";
            my @tmp;
            while( <$fhi> ) {
                chomp;
                next if /^Start/;
                my @ps = split(" ");
                my ($rS, $rE) = @ps[0..1];
                @tmp = @ps[2..$#ps] if @ps > 2;
                my ($st1, $st2) = @tmp[0..1];
                my $poss = join(" ", @tmp[2..$#tmp]);
                print $fh join("\t", $chr, $rS, $rE, $st1, $st2, $poss)."\n";
            }
            print "$round [$sS..$sE]: done\n";
        } else {
            print "$round [$sS..$sE]: failed\n";
        }
    }
}

sub pipe_fastphase {
    my ($dirW, $fwin) = rearrange(['dir', 'fwin'], @_);
    my $d01 = dir( $dirW, "01_in");
    my $d05 = dir( $dirW, '05');
    my $f06 = file($dirW, "06_log.txt");
    run_fastphase(-indir=>$d01, -fwin=>$fwin, -outdir=>$d05, -flog=>$f06);
}
sub run_fastphase {
    my ($dirI, $fWin, $dirO, $fLog) = rearrange(["indir", 'fwin', "outdir", "flog"], @_);
    die "inDir is not there\n" unless -d $dirI;
    system("mkdir -p $dirO") unless -d $dirO;
    my $t = readTable(-in=>$fWin, -header=>1);
    my $fOutH = new IO::File $fLog, "w";
    for my $i (0..$t->nofRow-1) {
        my ($round, $cntSnp) = map {$t->elm($i, $_)} qw/round cntSnp/;
        next unless $cntSnp > 1;
        next unless $round =~ /^5/ && $round ge "5_416_1";
        my $fi = file($dirI, "$round.txt");
        die "$fi is not there\n" unless -s $fi;
        my $cmd = qq/fastphase -KL2 -KU10 -H50 -i -o$round $fi/;
        my $t0 = [gettimeofday];
#    die $cmd."\n";
        runCmd($cmd, 1);
        system("mv $round* $dirO");
        my $t1 = [gettimeofday];
        printf $fOutH "%s\t%d\t%d\n", $round, $cntSnp, tv_interval($t0, $t1);
    }
}

sub pipe_haploview {
    my ($dirW, $fwin) = rearrange(['dir', 'fwin'], @_);
    my $d01 = dir( $dirW, "01_in");
    my $d05 = dir( $dirW, '05');
    my $d21 =  dir($dirW, "21_haploview_in" );
#  fastphase2haploview(-indir1=>$d20, -indir2=>$d01, -outdir=>$d21, -fwin=>$f02);
    my $d22 =  dir($dirW, "22_haploview_out");
#  run_haploview(-indir=>$d01, -outdir=>$d22, -fwin=>$f02);
    my $f23 = file($dirW, '23_tSNP.txt');
#  extract_tagSNP(-indir1=>$d01, -indir2=>$d22, -fout=>$f23, -fwin=>$f02);
}
sub run_haploview {
    my ($dirI, $dirO, $fWin) = rearrange(["indir", 'outdir', 'fwin'], @_);
    die("inDir is not there\n") unless -d $dirI;
    system("mkdir -p $dirO") unless -d $dirO;
    my $t = readTable(-in=>$fWin, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($round, $cntSnp) = map {$t->elm($i, $_)} qw/round cntSnp/;
        next unless $cntSnp > 0;
        my $fIn1  = file($dirI, "$round\_haploview.txt");
        my $fIn2  = file($dirI, "$round\_haploview_loc.txt");
        my $fOut = file($dirO, $round);
        print "$round\n";
        runHaploview(-in1=>$fIn1, -in2=>$fIn2, -out=>$fOut);
    }
}
sub runHaploview {
    my ($fi1, $fi2, $fo) = rearrange(["in1", "in2", "out"], @_);
    die("$fi1 is not there\n") unless -s $fi1;
    die("$fi2 is not there\n") unless -s $fi2;
    my $haploview = file($DIR_Bin, "Haploview.jar");
    my $fTmp = "/tmp/haploview.tmp";
    my $cmd = join(" ", "java", "-jar", $haploview, "-memory", 2048, "-nogui",
        "-pedfile", $fi1, "-info", $fi2, -out, $fTmp, "-hwcutoff", 0, "-maxDistance", 100,
        "-blockoutput", "GAB", "-compressedpng", "-pairwiseTagging");
    runCmd($cmd, 1);
    my $h = {
        'GABRIELblocks' => 'blocks.txt',
        'LD.PNG'        => 'LD.png',
        'TAGS'          => 'tags.txt',
        'TESTS'         => 'tests.txt'
    };
    for my $key (%$h) {
        my $fi = "$fo.$key";
        die("$fi is not there\n") unless -s $fi;
        move($fi, $fo."_".$h->{$key});
    }
}
sub extract_tagSNP {
    my ($dirI1, $dirI2, $fOut, $fWin) = rearrange(["indir1", "indir2", "fout", "fwin"], @_);
    my $t = readTable(-in=>$fWin, -header=>1);
    my $fOutH = new IO::File $fOut, "w";
    my $posPrev = 0;
    for my $i (0..$t->nofRow-1) {
        my ($round, $cntSnp) = map {$t->elm($i, $_)} qw/round cntSnp/;
        next unless $cntSnp > 0;
        my $fIn1 = file($dirI1, "$round\_haploview.txt");
        my $fIn2 = file($dirI2, "$round\_tests.txt");
        die("$fIn1 is not there\n") unless -s $fIn1;
        die("$fIn2 is not there\n") unless -s $fIn2;
        my $fInH1 = new IO::File $fIn1, "r";
        my $fInH2 = new IO::File $fIn2, "r";
        my $locHash = {};
        while( <$fInH1> ) {
            chomp;
            next if !$_;
            my @eleAry = split("\t", $_);
            $locHash->{$eleAry[0]} = $eleAry[1];
        }
        my (@tagSNP, @tagSNPUsed);
        while( <$fInH2> ) {
            chomp;
            next if !$_;
            my $locus = $_;
            $locus =~ s/\W//g;
            die("$locus not exist in $fInH1\n") unless exists $locHash->{$locus};
            push (@tagSNP, $locHash->{$locus});
            push (@tagSNPUsed, $locHash->{$locus}) if $locHash->{$locus} > $posPrev;
        }
        @tagSNPUsed = sort {$a<=>$b} @tagSNPUsed;
        print "input file ".sprintf("%03d", $round).": ".@tagSNPUsed." / ".@tagSNP." tagSNPs selected\n";
        print $fOutH join("\n", @tagSNPUsed, "\n");
        $posPrev = max(values(%$locHash));
    }
}

sub pipe_maxhap {
    my ($dirW, $fwin) = rearrange(['dir', 'fwin'], @_);
    system("mkdir -p $dirW") unless -d $dirW;
    my $d01 = dir( $dirW, "../02_simple_snp");
    my $f03 = file($dirW, "h26rho");
    my $d05 = dir( $dirW, '05');
    run_maxhap(-fwin=>$fwin, -indir=>$d01, -in2=>$f03, -outdir=>$d05);
}
sub run_maxhap {
    my ($dirI, $fcfg, $fWin, $dirO) = rearrange(['indir', 'in2', 'fwin', 'outdir'], @_);
    system("mkdir -p $dirO") unless -d $dirO;
    my $t = readTable(-in=>$fWin, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($round, $cntSnp) = map {$t->elm($i, $_)} qw/round cntSnp/;
        next unless $cntSnp > 1;
        next unless $round =~ /^[235]\_/;
        print "$round...\n";
        my $fi = file($dirI, "$round.txt");
        die "$fi is not there\n" unless -s $fi;
        my $fo = file($dirO, "$round.txt");
        my $cmd = qq/awk '{ if(NR==3) {print "anc", \$0} else {print \$0} }' $fi | exhap | maxhap 1 $fcfg 0.001 100 0.001 0 10 11 400 > $fo/;
        runCmd($cmd, 0);
    }
}

sub run_meme {
    my ($fIn, $fOut, $param) = rearrange(['in', 'out', 'param'], @_);
    unless($fOut) {
        my @tmp = reverse split(/\./, $fIn, 2);
        $fOut = $tmp[0];
    }
    my @opt = ("-oc $fOut");
    for my $key (keys %$param) {
        my $value = $param->{$key};
        if($key eq "flag") {
            my $aryRef = ref($value) eq "ARRAY" ? $value : [$value];
            for my $flag (@$aryRef) {
                push @opt, "-$flag";
            }
        } else {
            push @opt, "-$key $value";
        }
    }
    my $cmd = join(" ", "meme", $fIn, @opt);
    runCmd($cmd, 1);
}
sub run_mummer {
    my ($seqs, $fo) = rearrange(['seqs', 'out'], @_);
    my $dirE = "/soft/sle11/mummer/3.22";
    my ($f1, $f2) = @$seqs;
    my $f3 = "out.delta";
    runCmd("$dirE/nucmer --maxmatch -c 40 -p out $f1 $f2", 1);
    my $f4 = "out.coords";
    runCmd("$dirE/show-coords -r -c -l $f3 > $f4", 1);
    parse_mummer_coords($f4, $fo, 2);
    system("rm $f3 $f4");
}


1;
__END__
