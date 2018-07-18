#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use InitPath;
use Common; 
use Bio::Seq;
use Bio::SeqIO;
use Graph;
use Bio::SeqFeature::Generic;
use Readfile;
use Writefile;
use Annotate;
use Align;
use Qry;
use Getopt::Long;
use Time::HiRes qw/gettimeofday tv_interval/;
use Data::Dumper;
use Path::Class; 
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
sub readHapSample {
    my ($fi) = @_;
    my $fh = new IO::File $fi, "r";
    my $line = readline($fh);
    my @ps = split("\t", $line);
    my $nSample = @ps - 4;
    seek($fh, 0, 0);
    my (@markers, @gts);
    while(<$fh>) {
        chomp;
        next unless $_;
        my @ps = split "\t";
        die "not $nSample samples at:\n$_\n" unless @ps == $nSample + 4;
        push @markers, $ps[1];
        push @gts, [ @ps[4..$#ps] ];
    }
    my $nMarker = @gts;
    return ($nMarker, $nSample, \@gts, \@markers);
}
sub hapsample2Team {
#download_mdr_data();
    my $dirW = dir($DIR_In, "csci8980/01");
    my $fi1 = file($dirW, "case_genotypes_28409.dat");
    my $fi2 = file($dirW, "anticase_genotypes_28409.dat");
    my $fo1 = file($dirW, "team_geno.txt");
    my $fo2 = file($dirW, "team_pheno.txt");
    my ($nMarker1, $nSample1, $gts1, $markers1) = readHapSample($fi1);
    my ($nMarker2, $nSample2, $gts2, $markers2) = readHapSample($fi2);
    die "not equal markers: $nMarker1 != $nMarker2\n" unless $nMarker1 == $nMarker2; 
    printf "%4d markers, %4d samples (%4d cases + %4d controls)\n", $nMarker1, $nSample1+$nSample2, $nSample1, $nSample2;
    my $fho1 = new IO::File $fo1, "w";
    for my $i (0..$nMarker1-1) {
        my $row1 = $gts1->[$i];
        my $row2 = $gts2->[$i];
        print $fho1 join("", @$row1, @$row2)."\n";
    }
    my $fho2 = new IO::File $fo2, "w";
    print $fho2 join("", (1) x $nSample1, (2) x $nSample2)."\n";
}

my ($beg, $end) = (0, -1);
GetOptions('beg=i'=>\$beg, 'end=i'=>\$end);

my $dirW = dir($DIR_In, "csci8980/03_mdr_data");
#run_mdr_batch($dirW);
#sum_mdr($dirW);
sub download_mdr_data {
    my ($dirW) = @_;
    my @mafs = (0.2, 0.4);
    my @sizes = (200, 400, 800, 1600);
    my @hs = (0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4);
    my $url_pre = "http://discovery.dartmouth.edu/epistatic_data/";
    my $fn_format = "epistatic_data_%d_%s_%s.zip";
    for my $maf (@mafs) {
        for my $size (@sizes) {
            for my $h (@hs) {
                my $fn = sprintf $fn_format, $size, $maf, $h;
                my $url = $url_pre.$fn;
                system("wget -P $dirW $url");
                my $f = file($dirW, $fn);
                system("unzip $f -d $dirW");
                system("rm $f");
            }
        }
    }
}
sub parse_mdr_info {
#parse_mdr_info($dirW);
    my ($dir) = @_;
    my $fi = file($dir, "00_info.txt");
    my $fo = file($dir, "01_in.tbl");
    my @sizes = (200, 400, 800, 1600);
    my $t = Data::Table->new([], [qw/model maf h size/]);
    my $fhi = new IO::File $fi, "r";
    while(<$fhi>) {
        if(/^Model: (\d+)/) {
            my $model = $1;
            my @lines;
            for (1..5) {
                my $line = readline($fhi);
                chomp($line);
                push @lines, $line;
            }
            my $maf = $1 if $lines[0] =~ /^Minor Allele Frequency: ([\d\.]+)/;
            my $h = $1 if $lines[1] =~ /^Heritability: ([\d\.]+)/;
            my @pts = map {split " ", $_} @lines[2..4];
            my $pt = join(" ", @pts);
            
            for (@sizes) {
                $t->addRow([$model, $maf, $h, $_]);
            }
        }
    }
    my $fho = new IO::File $fo, "w";
    print $fho $t->tsv(1);
}
sub run_mdr_batch {
    my ($dir) = @_;
    my $fi = file($dir, "01_in.tbl");
    my $t = readTable(-in=>$fi, -header=>1);
    my $dirI = dir($dir, "01_mdr");

    my $rep = 100;
    my ($rowBeg, $rowEnd) = (276, 279);
    for my $i ($rowBeg..$rowEnd) {
        my ($model, $maf, $h, $size) = $t->row($i);
        
        my $fo = file($dir, "21_mdr", sprintf "%03d.tbl", $i);
        my $fh = new IO::File $fo, "w";
        print $fh join("\t", qw/model size rep time result/)."\n";

        for my $j (0..$rep-1) {
            my $fi = file($dirI, $size, $model, "$model.$size.$j.txt");
            die "$fi is not there\n" unless -s $fi;
            my ($rst, $time) = run_mdr($fi);
            my $rstStr = join(" ", map {join("-", @$_)} @$rst);
            print $fh join("\t", $model, $size, $j, $time, $rstStr)."\n";
            print "row[$i] model[$model] size[$size] rep[$j] time[$time]\n";
        }
    }
}
sub run_mdr {
    my ($fi) = @_;
    my $jar = $ENV{'src'}."/mdr-2.0_beta_8.4/mdr.jar";
    my $cmd = "java -jar $jar -cv=10 -min=2 -max=2 -parallel $fi";
    my $t1 = [gettimeofday];
    open(my $out, $cmd." 2>&1 |") or die "cannot run cmd: $cmd\n";
    my $tag = 0;
    my @rst;
    my $time;
    while(<$out>) {
        chomp;
        if(/^\#\#\# Top Models/) {
            $tag = 1;
        }
        if($tag == 1 && /^\s*$/) {
            $tag = 0;
        }
        if($tag == 1 && /^X/) {
            my @ps = split("\t", $_);
            push @rst, \@ps;
        }
    }
    my $t2 = [gettimeofday];
    my $tv = tv_interval($t1, $t2);
    return (\@rst, sprintf "%.02f", $tv);
}
sub sum_mdr {
    my ($dir) = @_;
    my $fw = file($dir, "01_in.tbl");
    my $tw = readTable(-in=>$fw, -header=>1);
    my $fo = file($dir, "21_mdr.tbl");
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/model maf h size time power1 power2/)."\n";

    for my $i (0..$tw->nofRow-1) {
        my ($model, $maf, $h, $size) = $tw->row($i);
        my $fi = file($dir, "21_mdr", sprintf "%03d.tbl", $i);
        my $t = readTable(-in=>$fi, -header=>1);
        my $time = sprintf "%.03f", sum($t->col("time")) / $t->nofRow;
        my ($cnt1, $cnt2) = (0, 0);
        for my $j (0..$t->nofRow-1) {
            my @ps = split(" ", $t->elm($j, "result"));
            my @ary;
            for my $str (@ps) {
                my ($str2, $score) = split("-", $str);
                my ($id1, $id2) = split(",", $str2);
                ($id1, $id2) = ($id2, $id1) if $id1 gt $id2;
                push @ary, [$id1, $id2, $score];
            }
            @ary = reverse sort {$a->[2] <=> $b->[2]} @ary;
            my $idx = first_index {$_->[0] eq "X0" && $_->[1] eq "X1"} @ary;
            $cnt1 ++ if $idx == 0;
            $cnt2 ++ if $idx >= 0 && $idx <= 2;
        }
        print $fh join("\t", $model, $maf, $h, $size, $time, $cnt1, $cnt2)."\n";
    }
}

#get_plink_input($dirW);
#run_plink_batch($dirW);
#sum_plink($dirW);
sub get_plink_input {
    my ($dir) = @_;
    my $fi = file($dir, "01_in.tbl");
    my $t = readTable(-in=>$fi, -header=>1);
    my $rep = 100;
    for my $i (0..$t->nofRow-1) {
        my ($model, $maf, $h, $size) = $t->row($i);
        my $dirI = dir($dir, "01_mdr", $size, $model);
        my $dirO = dir($dir, "02_plink", $size, $model);
        system("mkdir -p $dirO") unless -d $dirO;
        for my $j (0..$rep-1) {
            my $fi = file($dirI, "$model.$size.$j.txt");
            die "$fi is not there\n" unless -s $fi;
            my $fo = file($dirO, "$model.$size.$j.ped");
            convert_mdr_plink($fi, $fo);
        } 
    }
}
sub convert_mdr_plink {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $fh = new IO::File $fo, "w";
    for my $i(0..$t->nofRow-1) {
        my @ps = $t->row($i);
        my $id = sprintf "%03d", $i+1;
        my @gts = map {$_ == 0 ? "A A" : $_ == 1 ? "A T" : "T T"} @ps[0..$#ps-1];
        my $gt_str = join("  ", @gts);
        print $fh join("\t", $id, $id, -9, -9, -9, $ps[$#ps], $gt_str)."\n";
    }
}
sub run_plink_batch {
    my ($dir) = @_;
    my $fi = file($dir, "01_in.tbl");
    my $t = readTable(-in=>$fi, -header=>1);
    my $f_map = file($dir, "02_plink/1k.map");
    die "$f_map is not there\n" unless -s $f_map;
    my $rep = 100;
    for my $i ($beg..$end) {
        my ($model, $maf, $h, $size) = $t->row($i);
        my $fo = file($dir, "22_plink", sprintf "%03d.tbl", $i);
        my $fh = new IO::File $fo, "w";
        print $fh join("\t", qw/model size rep time result/)."\n";
        
        my $dirI = dir($dir, "02_plink", $size, $model);
        for my $j (0..$rep-1) {
            my $fi = file($dirI, "$model.$size.$j.ped");
            die "$fi is not there\n" unless -s $fi;
            my ($rst, $time) = run_plink($fi, $f_map, "$dir/22_plink/$model.$size.$j");
            my $rstStr = join(" ", map {join("-", @$_)} @$rst);
            print $fh join("\t", $model, $size, $j, $time, $rstStr)."\n";
            print "row[$i] model[$model] size[$size] rep[$j] time[$time]\n";
        } 
        close $fh;
    }
}
sub run_plink {
    my ($fi, $f_map, $pre) = @_;
    my $cmd = "plink --ped $fi --1 --allow-no-sex --map $f_map --map3 --fast-epistasis --out $pre";
    my $t1 = [gettimeofday];
    runCmd($cmd);
    my @rst;
    my $fh = new IO::File "$pre.epi.cc", "r";
    my $cnt = 0;
    while(<$fh>) {
        chomp;
        my @ps = split " ";
        next if ++$cnt == 1;
        push @rst, [ @ps[1,3,4,5] ];
        if($cnt > 300) {
            print "   > 300 pairs, remaining skipped\n";
            last;
        }
    }
    close $fh;
    system("rm $pre.epi.cc $pre.epi.cc.summary $pre.nof $pre.nosex $pre.log");
    my $t2 = [gettimeofday];
    my $tv = tv_interval($t1, $t2);
    return (\@rst, sprintf "%.02f", $tv);
} 
sub sum_plink {
    my ($dir) = @_;
    my $fw = file($dir, "01_in.tbl");
    my $tw = readTable(-in=>$fw, -header=>1);
    my $fo = file($dir, "22_plink.tbl");
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/model maf h size time power1 power2/)."\n";

    for my $i (0..$tw->nofRow-1) {
        my ($model, $maf, $h, $size) = $tw->row($i);
        my $fi = file($dir, "22_plink", sprintf "%03d.tbl", $i);
        next unless -s $fi;
        my $t = readTable(-in=>$fi, -header=>1);
        my $time = sprintf "%.03f", sum($t->col("time")) / $t->nofRow;
        my ($cnt1, $cnt2) = (0, 0);
        for my $j (0..$t->nofRow-1) {
            my @ps = split(" ", $t->elm($j, "result"));
            my @ary = map {[split("-", $_)]} @ps;
            @ary = reverse sort {$a->[2] <=> $b->[2]} @ary;
            my $idx = first_index {$_->[0] eq "X0" && $_->[1] eq "X1"} @ary;
            $cnt1 ++ if $idx == 0;
            $cnt2 ++ if $idx >= 0 && $idx <= 2;
        }
        print $fh join("\t", $model, $maf, $h, $size, $time, $cnt1, $cnt2)."\n";
    }
}

#get_beam_input($dirW);
#run_beam_batch($dirW);
sub get_beam_input {
    my ($dir) = @_;
    my $fi = file($dir, "01_in.tbl");
    my $t = readTable(-in=>$fi, -header=>1);
    my $rep = 100;
    for my $i (0..$t->nofRow-1) {
        my ($model, $maf, $h, $size) = $t->row($i);
        my $dirI = dir($dir, "01_mdr", $size, $model);
        my $dirO = dir($dir, "03_beam", $size, $model);
        system("mkdir -p $dirO") unless -d $dirO;
        for my $j (0..$rep-1) {
            my $fi = file($dirI, "$model.$size.$j.txt");
            die "$fi is not there\n" unless -s $fi;
            my $fo = file($dirO, "$model.$size.$j.txt");
            convert_mdr_beam($fi, $fo);
        } 
    }
}
sub convert_mdr_beam {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $fh = new IO::File $fo, "w";
    my $n_snp = $t->nofCol - 1;
    my $n_ind = $t->nofRow;
    my @pts = map {$t->elm($_, $n_snp)} (0..$n_ind-1);
    print $fh join("\t", qw/ID Chr Pos/, join(" ", @pts))."\n";
    my @snp_ids = $t->header;
    for my $i (0..$n_snp-1) {
        my $snp_id = $snp_ids[$i];
        my $pos = $i * 10000 + 1;
        my @gts = map {$t->elm($_, $snp_id)} (0..$n_ind-1);
        my $gt_str = join(" ", @gts);
        print $fh join("\t", $snp_id, "chr1", $pos, $gt_str)."\n";
    }
}

#get_team_input($dirW);
#run_team_batch($dirW);
#sum_team($dirW);
sub get_team_input {
    my ($dir) = @_;
    my $fi = file($dir, "01_in.tbl");
    my $t = readTable(-in=>$fi, -header=>1);
    my $rep = 100;
    for my $i (240..$t->nofRow-1) {
        my ($model, $maf, $h, $size) = $t->row($i);
        my $dirI = dir($dir, "01_mdr", $size, $model);
        my $dirO = dir($dir, "04_team", $size, $model);
        system("mkdir -p $dirO") unless -d $dirO;
        for my $j (0..$rep-1) {
            my $fi = file($dirI, "$model.$size.$j.txt");
            die "$fi is not there\n" unless -s $fi;
            my $fo = file($dirO, "$model.$size.$j");
            convert_mdr_team($fi, $fo);
        } 
    }
}
sub convert_mdr_team {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my ($fo1, $fo2) = map {"$fo\_$_.txt"} qw/gt pt/;
    my $fh1 = new IO::File $fo1, "w";
    my $fh2 = new IO::File $fo2, "w";
    my $n_snp = $t->nofCol - 1;
    my $n_ind = $t->nofRow;

    my @pts = map {$t->elm($_, $n_snp)} (0..$n_ind-1);
    print $fh2 join("", @pts)."\n";

    my @snp_ids = $t->header;
    for my $i (0..$n_snp-1) {
        my $snp_id = $snp_ids[$i];
        my @gts = map {$t->elm($_, $snp_id)} (0..$n_ind-1);
        print $fh1 join("", @gts)."\n";
    }
}
sub run_team_batch {
    my ($dir) = @_;
    my $fi = file($dir, "01_in.tbl");
    my $t = readTable(-in=>$fi, -header=>1);
    my $rep = 100;
    for my $i ($beg..$end) {
        my ($model, $maf, $h, $size) = $t->row($i);
        my $fo = file($dir, "24_team", sprintf "%03d.tbl", $i);
        my $fh = new IO::File $fo, "w";
        print $fh join("\t", qw/model size rep time result/)."\n";
        
        my $dirI = dir($dir, "04_team", $size, $model);
        for my $j (0..$rep-1) {
            my ($fi1, $fi2) = map {file($dirI, "$model.$size.$j\_$_.txt")} qw/gt pt/;
            die "$fi1 is not there\n" unless -s $fi1;
            die "$fi2 is not there\n" unless -s $fi2;
            my ($rst, $time) = run_team($fi1, $fi2, $size, 1000, 100, 0.2);
            my $rstStr = join(" ", map {join("-", @$_)} @$rst);
            print $fh join("\t", $model, $size, $j, $time, $rstStr)."\n";
            print "row[$i] model[$model] size[$size] rep[$j] time[$time]\n";
        } 
        close $fh;
    }
}
sub run_team {
    my ($fi1, $fi2, $n_ind, $n_snp, $n_per, $fde) = @_;
    my $exe = $ENV{'src'}."/TEAM-0.0.4_linux/TEAM/team.sh";
    my $cmd = "$exe $fi1 $fi2 $n_ind $n_snp $n_per $fde";
    my $t1 = [gettimeofday];
    open(my $out, $cmd." 2>&1 |") or die "cannot run cmd: $cmd\n";
    my $tag = 0;
    my @rst;
    my $time;
    while(<$out>) {
        chomp;
        if(/^\=+Significant SNP\-pairs/) {
            $tag = 1;
        }
        if($tag == 1 && /^\=+$/) {
            $tag = 0;
        }
        if($tag == 1 && /^SNP\-pair \((\d+),(\d+)\)\:([\d\.]+)/) {
            push @rst, [$1, $2, $3];
        }
    }
    my $t2 = [gettimeofday];
    my $tv = tv_interval($t1, $t2);
    return (\@rst, sprintf "%.02f", $tv);
}
sub sum_team {
    my ($dir) = @_;
    my $fw = file($dir, "01_in.tbl");
    my $tw = readTable(-in=>$fw, -header=>1);
    my $fo = file($dir, "24_team.tbl");
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/model maf h size time power1 power2/)."\n";

    for my $i (0..$tw->nofRow-1) {
        my ($model, $maf, $h, $size) = $tw->row($i);
        my $fi = file($dir, "24_team", sprintf "%03d.tbl", $i);
        next unless -s $fi;
        next if $size == 200;
        my ($time, $cnt1, $cnt2) = parse_team_output($fi);
        print join("\t", $model, $maf, $h, $size, $time, $cnt1, $cnt2)."\n";
        print $fh join("\t", $model, $maf, $h, $size, $time, $cnt1, $cnt2)."\n";
    }
}
sub parse_team_output {
    my ($fi) = @_;
    my $fh = new IO::File $fi, "r";
    my @times;
    my ($cnt1, $cnt2) = (0, 0);
    while(<$fh>) {
        chomp;
        my @ps = split "\t";
        next if $ps[0] eq "model";
        push @times, $ps[3];
        next unless $ps[4];
        my @hits = split(" ", $ps[4]);
        my @ary = map {[split "-"]} @hits;
        @ary = reverse sort {$a->[2] <=> $b->[2]} @ary;
        my $idx = first_index {$_->[0] == 0 && $_->[1] == 1} @ary;
        $cnt1 ++ if $idx == 0;
        $cnt2 ++ if $idx >= 0 && $idx <= 2;
    }
    my $time = sprintf "%.03f", sum(@times) / @times;
    $cnt1 = sprintf "%d", $cnt1 * 100/@times;
    $cnt2 = sprintf "%d", $cnt2 * 100/@times;
    return ($time, $cnt1, $cnt2);
}

#get_boost_input($dirW);
#run_boost_batch($dirW);
#sum_boost($dirW);
sub get_boost_input {
    my ($dir) = @_;
    my $fi = file($dir, "01_in.tbl");
    my $t = readTable(-in=>$fi, -header=>1);
    my $rep = 100;
    for my $i (0..$t->nofRow-1) {
        my ($model, $maf, $h, $size) = $t->row($i);
        my $dirI = dir($dir, "01_mdr", $size, $model);
        my $dirO = dir($dir, "05_boost", $size, $model);
        system("mkdir -p $dirO") unless -d $dirO;
        for my $j (0..$rep-1) {
            my $fi = file($dirI, "$model.$size.$j.txt");
            die "$fi is not there\n" unless -s $fi;
            my $fo = file($dirO, "$model.$size.$j.txt");
            convert_mdr_boost($fi, $fo);
        } 
    }
}
sub convert_mdr_boost {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $fh = new IO::File $fo, "w";
    my $n_snp = $t->nofCol - 1;
    my $n_ind = $t->nofRow;
    for my $i (0..$n_ind-1) {
        my @ps = $t->row($i);
        print $fh join("\t", @ps[$#ps, 0..$#ps-1])."\n";
    }
    close $fh;
}
sub run_boost_batch {
    my ($dir) = @_;
    my $fi = file($dir, "01_in.tbl");
    my $t = readTable(-in=>$fi, -header=>1);
    my $rep = 100;
    my ($rowBeg, $rowEnd) = (0, 279);
    for my $i ($rowBeg..$rowEnd) {
        my ($model, $maf, $h, $size) = $t->row($i);
        my $fo = file($dir, "25_boost", sprintf "%03d.tbl", $i);
        my $fh = new IO::File $fo, "w";
        print $fh join("\t", qw/model size rep time result/)."\n";
        
        my $dirI = dir($dir, "05_boost", $size, $model);
        for my $j (0..$rep-1) {
            my $fi = file($dirI, "$model.$size.$j.txt");
            die "$fi is not there\n" unless -s $fi;
            my ($rst, $time) = run_boost($fi, $dir);
            my $rstStr = join(" ", map {join(",", @$_)} @$rst);
            print $fh join("\t", $model, $size, $j, $time, $rstStr)."\n";
            print "row[$i] model[$model] size[$size] rep[$j] time[$time]\n";
        } 
        close $fh;
    }
}
sub run_boost {
    my ($fi, $dir) = @_;
    my $dirW = dir($dir, "05_boost");
    chdir($dirW);
    my $ftmp = file($dirW, "filenamelist.txt");
    system("echo $fi > $ftmp");
    my $t1 = [gettimeofday];
    runCmd("BOOST");
    my $t2 = [gettimeofday];
    my $tv = tv_interval($t1, $t2);

    my @rst, 
    my $fo = file($dirW, "InteractionRecords.txt");
    my $fh = new IO::File $fo, "r";
    while(<$fh>) {
        chomp;
        my @ps = split(" ", $_);
        push @rst, [ @ps[1..2], map {sprintf "%.02f", $_} @ps[3..6] ];
    }
    return (\@rst, sprintf "%.02f", $tv);
}
sub sum_boost {
    my ($dir) = @_;
    my $fw = file($dir, "01_in.tbl");
    my $tw = readTable(-in=>$fw, -header=>1);
    my $fo = file($dir, "25_boost.tbl");
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/model maf h size time power1 power2/)."\n";

    for my $i (0..$tw->nofRow-1) {
        my ($model, $maf, $h, $size) = $tw->row($i);
        my $fi = file($dir, "25_boost", sprintf "%03d.tbl", $i);
        my $t = readTable(-in=>$fi, -header=>1);
        my $time = sprintf "%.03f", sum($t->col("time")) / $t->nofRow;
        my ($cnt1, $cnt2) = (0, 0);
        for my $j (0..$t->nofRow-1) {
            my @ps = split(" ", $t->elm($j, "result"));
            next unless @ps;
            my @ary = map {[split ","]} @ps;
            @ary = reverse sort {$a->[4] <=> $b->[4]} @ary;
            my $idx = first_index {$_->[0] == 0 && $_->[1] == 1} @ary;
            $cnt1 ++ if $idx == 0;
            $cnt2 ++ if $idx >= 0 && $idx <= 2;
        }
        print $fh join("\t", $model, $maf, $h, $size, $time, $cnt1, $cnt2)."\n";
    }
}


sub genomeSim2Mdr {
    my ($dirW) = @_;
    my ($dirI, $dirO) = map {dir($dirW, $_)} qw/data data_mdr/;
    system("mkdir -p $dirO") unless -d $dirO;
    for my $i (1..10) {
        my $fn = sprintf "my.0.1.%d-cc1000.mdr", $i;
        my $fi = file($dirI, $fn);
        my $fo = file($dirO, sprintf("%02d.txt", $i));
        my $fhi = new IO::File $fi, "r";
        my $fho = new IO::File $fo, "w";
        my $l = 0;
        while(<$fhi>) {
            chomp;
            my @ps = split " ";
            if($l++ == 0) {
                my @labels = map {sprintf "X%02d", $_} (1..$#ps);
                print $fho join("\t", @labels, $ps[0])."\n";
            }
            print $fho join("\t", @ps[1..$#ps], $ps[0])."\n";
        }
    }
}
$dirW = dir($DIR_In, "csci8980/genomeSIMLA/01");
#genomeSim2Mdr($dirW);


