#!/usr/bin/perl
use strict;
use InitPath;
use Common;
use WindowStat;
use Medicago;
use Data::Dumper;
use Path::Class; 

my $refDb = "mt_35";
my $opt = "acc288";
my $ids = get_acc_ids($opt);

my $dirW = dir($DIR_Repo, $refDb, "30_vnt_$opt", "21_coverage");
my $chrs = [map { "chr".$_ } (1..8)];
my $d01 = dir($dirW, "01_raw");
my $d02 = dir($dirW, "02");
#processCov($d01, $d02, $ids);

my $d11 = dir($dirW, "11_chr_pos");
#createChrPosFile($d11, $chrs, $refDb);
#mergeCovFiles($dirW, $ids, $chrs, 1, $refDb);
#mergeCovFiles($dirW, $ids, $chrs, 2, $refDb);
#writeCovInfo($dirW, $ids, $chrs, $refDb);
#print "cov_window -i $dirW/21_windows.tbl -o $dirW/22_cov_info.tbl -t acc288 -c 2\n";

sub processCov {
    my ($d01, $d02, $ids) = @_;
    my $dirI = dir($d01, "coverage_unique");
    my $dirO = dir($d02, "coverage_unique");
    system("mkdir -p $dirO") unless -d $dirO;

    my %chrH = map { "chr".$_ => "Mt3.5.1Chr".$_} (1..2, 4..8);
    $chrH{"chr3"} = "Mt3.5Chr3";
    for my $id (@$ids) {
#    next unless $id eq "HM020";
        my $f_gz = file($dirI, "$id.tgz");
        my $lines = runCmd2("tar -tf $f_gz | head -n 1");
        my @ps = split("/", $lines->[0]);
        my $pre = $ps[0];

        my $fo_str = join(" ", map {"$pre/".$_} values(%chrH));
        runCmd("tar -C $dirO -xzf $f_gz $fo_str", 1);
        if($id ne $pre) {
            runCmd("rm -rf $dirO/$id") if -d "$dirO/$id";
            runCmd("mv $dirO/$pre $dirO/$id");
        } 
        print "\trenaming chromosome files\n";
        for my $chr (keys(%chrH)) {
            system("mv $dirO/$id/".$chrH{$chr}." $dirO/$id/$chr");
        }
    }
}
sub createChrPosFile {
    my ($dir, $chrs, $refDb) = @_;
    system("mkdir -p $dir") unless -d $dir;
    print "creating chr position files:\n";
    for my $chr (@$chrs) {
        print "\tworking on $chr\n";
        my $f = file($dir, "$chr.tbl");
        my $fh = new IO::File $f, "w";
        for my $i (1..getSeqLen($chr, $refDb)) {
            print $fh join("\t", $chr, $i, $i)."\n";
        }
    }
}
sub mergeCovFiles {
    my ($dirW, $ids, $chrs, $opt, $refDb) = @_;
    my %chrs;
    if($refDb eq "mt_30") {
        %chrs = map { "chr".$_ => "MtChr".$_ } (1..8);
        $chrs{"chrCp"} = "chloroplast_AC093544.8";
    } elsif($refDb eq "mt_35") {
        %chrs = map { "chr".$_ => "Mt3.5.1Chr".$_ } (1..2,4..8);
        $chrs{"chr3"} = "Mt3.5Chr3";
#    $chrs{"chrCp"} = "chl_Mt";
    }

    my $d01 = dir($dirW, "01_raw");
    my $d11 = dir($dirW, "11_chr_pos");
    my $d12 = dir($dirW, "12_merged");
    my $d99 = dir($dirW, "99_tmp");
    my $t = readTable(-in=>file($d01, "pre.tbl"), -header=>1);
    my $h_pre_gz = { map {$t->elm($_, "id") => $t->elm($_, "pre")} (0..$t->nofRow-1) };
    my @f_chrs;
    my $pre_raw = $opt == 1 ? "coverage" : "coverage_unique";
    for my $chr (@$chrs) {
        print "\tworking on $chr:\n";
        my @fis = (file($d11, "$chr.tbl"));
        my $f_chr_name = $chrs{$chr};
        for my $id (@$ids) {  
            my $f_tgz = file($d01, $pre_raw, "$id.tgz");
            die "cannot open $f_tgz\n" unless -s $f_tgz;
            my $pre_gz = $h_pre_gz->{$id};
            runCmd("tar -C $d99 -xzf $f_tgz $pre_gz/$f_chr_name", 1); 
            my $f_cov = file($d99, $pre_gz, $f_chr_name);
            push @fis, $f_cov;
        }
        printf "merging coverage files for %s from %d accessions\n", $chr, scalar(@$ids);
        my $f_tbl = file($d12, "$chr.tbl.gz");
        runCmd("paste ".join(" ", @fis)." | gzip > $f_tbl", 1);
        system("rm -rf $d99/*");
        push @f_chrs, $f_tbl;
    }

    my $fo = file($dirW, "cov$opt.tbl.gz");
    runCmd("zcat ".join(" ", @f_chrs)." | bgzip > $fo", 1);
    runCmd("tabix -f -s1 -b2 -e3 $fo", 1);
    system("rm ".join(" ", @f_chrs));
}
sub writeCovInfo {
    my ($dir, $ids, $chrs, $refDb) = @_;
    my $fo  = file($dir, "21_windows.tbl");
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id chr beg end/)."\n";
    my $i = 0;
    for my $chr (@$chrs) {
        my $locAry = getWindows(-chr=>$chr, -winsize=>10000, -winstep=>10000, -db=>$refDb);
        for (@$locAry) {
            my ($wbeg, $wend) = @$_;
            print $fh join("\t", ++$i, $chr, $wbeg, $wend)."\n";
        }
    }
}



