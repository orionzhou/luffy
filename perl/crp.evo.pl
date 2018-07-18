#!/usr/bin/perl -w
use strict;
use lib $ENV{"SCRIPT_HOME_PERL"};
use File::Path qw/make_path remove_tree/;
use GD;
use GD::Image;
use GD::Polygon;
use Bio::Graphics;
use Bio::SeqFeature::Generic;
use Common;
use Seq;
use Align;
use BioFeature;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $org = "Athaliana";
my $dir_i = "/project/youngn/zhoup/Data/misc3/spada.crp/$org";
my $f_ref = "$dir_i/01_genome/01_refseq.fa";
my $f_gtb = "$dir_i/01_genome/61_gene.gtb";
my $f_pre = "$dir_i/31_model_SPADA/61_final.gtb";

my $dir = "/project/youngn/zhoup/Data/misc2/crp.evo/$org";
my $f01 = "$dir/01.gtb";
my $f_fam = "$dir/../../crp.annotation/family_info.tbl";
#add_fam_info($f_pre, $f_fam, $f01);
my $f02 = "$dir/02_seq.tbl";
#prepare_seq($f01, $f_ref, $f02);
my $f05 = "$dir/05_pairs.tbl";
#screen_pairs_raw($f02, $f05);
my $f07 = "$dir/07_yass.tbl";
#screen_pairs_yass($f02, $f05, $f07);
my $f09 = "$dir/09_yass_filtered.tbl";
#filter_yass($f07, $f09);
my $f11 = "$dir/11.tbl";
#prepare_plot($f02, $f09, $f11);
my $d21 = "$dir/21_figs";
#comparative_plot($f11, $d21, $f_gtb);

sub add_fam_info {
    my ($fgi, $ff, $fgo) = @_;
    my $tgi = readTable(-in=>$fgi, -header=>1);
    my $tf = readTable(-in=>$ff, -header=>1);
    my $h = { map {$tf->elm($_, "id") => $tf->elm($_, "cat")} (0..$tf->nofRow-1) };

    for my $i (0..$tgi->nofRow-1) {
        my $fam = $tgi->elm($i, "cat3");
        my $cat = $h->{$fam};
        $tgi->setElm($i, "note", $cat);
    }

    open(FH, ">$fgo") or die "cannot write to $fgo\n";
    print FH $tgi->tsv(1);
    close FH;
}
sub prepare_seq {
    my ($f_pre, $f_ref, $fo) = @_;
    my $t = readTable(-in=>$f_pre, -header=>1);

    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO join("\t", qw/id fam cat chr beg end strand locC seq_pro seq_cds seq_ext/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $cat, $chr, $begM, $endM, $str, $locCS, $seq_pro) = 
          map {$t->elm($i, $_)} qw/id cat3 note chr beg end strand locC seq/;
        my $locC = locStr2Ary($locCS);
        my $seq_cds = seqRet($locC, $chr, $str, $f_ref);

        my $beg = max(1, $begM-1000);
        my $end = min($endM+1000, seqLen($chr, $f_ref));
        my $seq_ext = seqRet([[$beg, $end]], $chr, $str, $f_ref);
        $locC = [ sort {$a->[0] <=> $b->[0]} @$locC ];
        $locC = [ reverse @$locC ] if $str eq "-";
        my $locCR = $str eq "-" ?
          [ map {[$end-$_->[1]+1, $end-$_->[0]+1]} @$locC ] : 
          [ map {[$_->[0]-$beg+1, $_->[1]-$beg+1]} @$locC ];
        print FHO join("\t", $id, $fam, $cat, $chr, $beg, $end, $str, locAry2Str($locCR),
          $seq_pro, $seq_cds, $seq_ext)."\n";
    }
    close FHO;
}
sub screen_pairs_raw {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO join("\t", qw/id1 id2 identity/)."\n";

    my $h = group($t->colRef("fam"));
    for my $fam (sort(keys(%$h))) {
        my ($idx, $cnt) = @{$h->{$fam}};
        for my $i ($idx..$idx+$cnt-1) {
            for my $j ($i+1..$idx+$cnt-1) {
                my ($id1, $chr1, $beg1, $end1, $seq1) = 
                  map {$t->elm($i, $_)} qw/id chr beg end seq_cds/;
                my ($id2, $chr2, $beg2, $end2, $seq2) = 
                  map {$t->elm($j, $_)} qw/id chr beg end seq_cds/;
#                next if ($chr1 eq $chr2 && abs($beg1-$beg2) < 100_000);
                my $seqO1 = Bio::Seq->new(-id=>$id1, -seq=>$seq1);
                my $seqO2 = Bio::Seq->new(-id=>$id2, -seq=>$seq2);
                my ($lenA, $lenM, $lenI, $lenG, $len1, $len2) = run_pw_aln($seqO1, $seqO2, 2);
#                die Dumper($lenA, $lenM, $lenI, $lenG) if $id1 eq "h0003.01" && $id2 eq "h0059.01";
                my $idty = sprintf "%.03f", $lenM/$lenA;
                next if $idty < 0.7;
                print FHO join("\t", $id1, $id2, $idty)."\n";
            }
        }
    }
    close FHO;
}
sub run_yass {
    my ($seq1, $seq2) = @_;
    my $f_bin = "/project/youngn/zhoup/Source/yass-1.14/src/yass";
    my $fi1 = "yass_seq1_".int(rand(1000)).".fa";
    my $fi2 = "yass_seq2_".int(rand(1000)).".fa";
    my $fo = "yass_".int(rand(1000)).".txt";
    writeFile($fi1, ">seq1", $seq1);
    writeFile($fi2, ">seq2", $seq2);
    my $cmd = "$f_bin $fi1 $fi2 -d 2 -o $fo";
    runCmd($cmd, 0);

    my @ary;
    open(FH, "<$fo") or die "cannot open $fo for reading\n";
    while(<FH>) {
        chomp;
        next if /^\#/;
        my @ps = split /\s+/;
        next if $ps[0] ne "seq1" || $ps[1] ne "seq2";
        die "yass output not 12 fields:\n$_\n$cmd\n" if @ps != 12;
        push @ary, [@ps[2..11]];
    }
    close FH;

    system("rm $fi1 $fi2 $fo");
    return \@ary;
}
sub screen_pairs_yass {
    my ($f01, $f02, $f03) = @_;
    
    my $hs;
    my $t1 = readTable(-in=>$f01, -header=>1);
    for my $i (0..$t1->nofRow-1) {
        my ($id, $locS, $seq_ext) = map {$t1->elm($i, $_)} qw/id locC seq_ext/;
        $hs->{$id} = [ locStr2Ary($locS), $seq_ext ];
    }

    open(FHO, ">$f03") or die "cannot open $f03 for writing\n";
    print FHO join("\t", qw/id1 id2 locC1 locC2 match/)."\n";

    my $t = readTable(-in=>$f02, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($id1, $id2) = $t->row($i);
        my ($loc1, $seq1) = @{$hs->{$id1}};
        my ($loc2, $seq2) = @{$hs->{$id2}};

        my $ary = run_yass($seq1, $seq2);
        my @aryStr;
        for (@$ary) {
            my ($pct_idty, $len_aln, $mm, $gap, $b1, $e1, $b2, $e2, $e, $score) = @$_;
            my $str = ($e2 < $b2) ? "-" : "+";
            ($b2, $e2) = ($e2, $b2) if $str eq "-";
            push @aryStr, join("_", $b1, $e1, $str, $b2, $e2, $pct_idty/100, $e, $len_aln);
        }
        my $str = join(",", @aryStr);
        print FHO join("\t", $id1, $id2, locAry2Str($loc1), locAry2Str($loc2), $str)."\n";
    }
    close FHO;
}
sub filter_yass {
    my ($f03, $f04) = @_;
    my $t = readTable(-in=>$f03, -header=>1);

    open(FH, ">$f04") or die "cannot open $f04 for writing\n";
    print FH join("\t", qw/id id1 id2 locC1 locC2 match/)."\n";
    my $id = 0;
    for my $i (0..$t->nofRow-1) {
        my ($id1, $id2, $locS1, $locS2, $matchStr) = $t->row($i);
        my $loc1 = locStr2Ary($locS1);
        my $loc2 = locStr2Ary($locS2);

        my $flag = 0;
        my @msf;
        my @ms = map {[split("_", $_)]} split(",", $matchStr);
        for (@ms) {
            my ($b1, $e1, $str, $b2, $e2, $pct, $e, $len) = @$_;
            next if $len < 50 || $e > 0.1 || $str eq "-";
            my ($locO1, $lenO1) = posOvlp($loc1, [[$b1, $e1]]);
            my ($locO2, $lenO2) = posOvlp($loc2, [[$b2, $e2]]);
            $flag = 1 if $lenO1 > 100 && $lenO2 > 100;
            $pct = sprintf "%.02f", $pct;
            push @msf, [$b1, $e1, $str, $b2, $e2, $pct, $e, $len] if $e <= 0.1 && $len >= 50;
        }
        next if $flag == 0;

        $matchStr = join(",", map {join("_", @$_)} @msf);
        print FH join("\t", ++$id, $id1, $id2, $locS1, $locS2, $matchStr)."\n";
    }
    close FH;
}
sub prepare_plot {
    my ($fs, $fy, $fo) = @_;
    my $ts = readTable(-in=>$fs, -header=>1);
    my $hs = { map {$ts->elm($_, "id") => $ts->rowRef($_)} (0..$ts->nofRow-1) };
    my $t = readTable(-in=>$fy, -header=>1);

    open(FH, ">$fo") or die "cannot open $fo for writing\n";
    print FH join("\t", qw/id id1 chr1 beg1 end1 str1 locC1 id2 chr2 beg2 end2 str2 locC2 match/)."\n";
    my $id = 0;
    for my $i (0..$t->nofRow-1) {
        my ($id, $id1, $id2, $locS1, $locS2, $matchStr) = $t->row($i);
        my ($chr1, $beg1, $end1, $str1) = @{$hs->{$id1}}[3..6];
        my ($chr2, $beg2, $end2, $str2) = @{$hs->{$id2}}[3..6];
        my $locCR1 = locStr2Ary($locS1);
        my $locCR2 = locStr2Ary($locS2);
        my $locC1 = $str1 eq "-" ? [ map {[$end1-$_->[1]+1, $end1-$_->[0]+1]} reverse(@$locCR1) ] : [ map {[$beg1+$_->[0]-1, $beg1+$_->[1]-1]} @$locCR1 ];
        my $locC2 = $str2 eq "-" ? [ map {[$end2-$_->[1]+1, $end2-$_->[0]+1]} reverse(@$locCR2) ] : [ map {[$beg2+$_->[0]-1, $beg2+$_->[1]-1]} @$locCR2 ];

        my @ps;
        for (split(",", $matchStr)) {
            my ($bl1, $el1, $strl, $bl2, $el2, $pct, $e, $len) = split("_", $_);
            my $strm1 = $str1 eq "+"   ? "+" : "-";
            my $strm2 = $str2 eq $strl ? "+" : "-";
            my ($bg1, $eg1) = $strm1 eq "+" ? ($beg1+$bl1-1, $beg1+$el1-1) : ($end1-$el1+1, $end1-$bl1+1);
            my ($bg2, $eg2) = $strm2 eq "+" ? ($beg2+$bl2-1, $beg2+$el2-1) : ($end2-$el2+1, $end2-$bl2+1);
            push @ps, join("_", $bg1, $eg1, $strm1, $bg2, $eg2, $strm2, $pct, $e, $len);
        }
        print FH join("\t", $id, $id1, $chr1, $beg1, $end1, $str1, locAry2Str($locC1), 
            $id2, $chr2, $beg2, $end2, $str2, locAry2Str($locC2), join(",", @ps))."\n";
    }
    close FH;
}
sub comparative_plot {
    my ($fi, $do, $f_gtb) = @_;
    make_path($do) unless -d $do;
    remove_tree($do, {keep_root=>1});
    my $t = readTable(-in=>$fi, -header=>1);
    my $tg = readTable(-in=>$f_gtb, -header=>1);
  
    for my $i (0..$t->nofRow-1) {
        my ($id, $id1, $chr1, $beg1, $end1, $str1, $locS1, $id2, $chr2, $beg2, $end2, $str2, $locS2, $matchS) = $t->row($i);
        my ($loc1, $loc2) = (locStr2Ary($locS1), locStr2Ary($locS2));
        my @ms = map {[split("_", $_)]} split(",", $matchS);
        next if $i > 10;

        my $p1 = Bio::Graphics::Panel->new( -start=>$beg1, -end=>$end1, -width=>800, -key_style=>'between', -grid=>1 );
        my $p2 = Bio::Graphics::Panel->new( -start=>$beg2, -end=>$end2, -width=>800, -key_style=>'between', -grid=>1 );
        my $ruler1 = Bio::SeqFeature::Generic->new( -display_name=>$chr1, -start=>$beg1, -end=>$end1 );
        my $ruler2 = Bio::SeqFeature::Generic->new( -display_name=>$chr2, -start=>$beg2, -end=>$end2 );
        
        my $fe1 = loc2Feature(-id=>$id1, -locC=>$loc1, -strand=>$str1);
        my $fe2 = loc2Feature(-id=>$id2, -locC=>$loc2, -strand=>$str2);

        my $fes_gene_1 = gtb2Features($tg, $chr1, $beg1, $end1);
        my $fes_gene_2 = gtb2Features($tg, $chr2, $beg2, $end2);
        
        my @fes_aln_1;
        my @fes_aln_2;
        my @ps;
        for my $k (1..@ms) {
            my ($b1, $e1, $str1, $b2, $e2, $str2, $pct, $e, $len) = @{$ms[$k-1]};
            my $fe_aln_1 = Bio::SeqFeature::Generic->new( -display_name=>"$id#$k",
                -start=>$b1, -end=>$e1, -strand=>$str1, -source=>$e, -score=>$pct );
            my $fe_aln_2 = Bio::SeqFeature::Generic->new( -display_name=>"$id#$k",
                -start=>$b2, -end=>$e2, -strand=>$str2, -source=>$e, -score=>$pct );
            push @fes_aln_1, $fe_aln_1;
            push @fes_aln_2, $fe_aln_2;
            my ($pb1, $pe1) = $p1->location2pixel($b1, $e1);
            my ($pb2, $pe2) = $p2->location2pixel($b2, $e2);
            push @ps, [$pb1, $pe1, $str1, $pb2, $pe2, $str2];
        }
 
        $p1->add_track($fes_gene_1, -glyph=>'processed_transcript', 
            -bgcolor=>'orange', -fgcolor=>'black', -connector=>'solid',
            -implied_utrs=> 1, -label=>1, -description => 0, -key=>'Genome Annotation'); 
        $p1->add_track($fe1, -glyph=>'processed_transcript',
            -bgcolor=>'skyblue', -fgcolor=>'slateblue', -connector=>'solid',
            -implied_utrs=> 1, -label=>1, -description => 0, -key=>'SPADA prediction'); 
        $p1->add_track($ruler1, -glyph=>'arrow', -double=>1, -tick=>2, -label=>1);
        $p1->add_track( \@fes_aln_1, -glyph=>'arrow', -fgcolor=>'purple',
            -label=>sub {my $fe = shift; 
                return $fe->display_name."[".$fe->source_tag."][".$fe->score."]";} );
        $p2->add_track( \@fes_aln_2, -glyph=>'arrow', -fgcolor=>'purple',
            -label=>0);
        $p2->add_track($ruler2, -glyph=>'arrow', -double=>1, -tick=>2, -label=>1);
        $p2->add_track($fe2, -glyph=>'processed_transcript',
            -bgcolor=>'skyblue', -fgcolor=>'slateblue', -connector=>'solid',
            -implied_utrs=> 1, -label=>1, -description => 0, -key=>'SPADA prediction'); 
        $p2->add_track($fes_gene_2, -glyph=>'processed_transcript',
            -bgcolor=>'orange', -fgcolor=>'black', -connector=>'solid',
            -implied_utrs=> 1, -label=>1, -description => 0, -key=>'Genome Annotation'); 
        
        my $gd1 = $p1->gd();
        my $gd2 = $p2->gd();
        my $wd = max($gd1->width(), $gd2->width());
        my $ht = $gd1->height() + 50 + $gd2->height();
        my $y1 = $gd1->height();
        my $y2 = $y1 + 50;

        my $gd = GD::Image->new($wd, $ht);
        my $co_wh = $gd->colorAllocate(255,255,255);
        my $co_bl = $gd->colorAllocate(  0,  0,  0);
        my $co_pe = $gd->colorAllocate(255,218,185);
        $gd->copy($gd1, 0,  0,0,0,$gd1->width(), $gd1->height());
        $gd->copy($gd2, 0,$y2,0,0,$gd2->width(), $gd2->height());
        for (@ps) {
            my ($pb1, $pe1, $str1, $pb2, $pe2, $str2) = @$_;
            my $pl = GD::Polygon->new();
            $pl->addPt($pb1, $y1);
            $pl->addPt($pe1, $y1);
            if($str1 eq $str2) {
                $pl->addPt($pe2, $y2);
                $pl->addPt($pb2, $y2);
            } else {
                $pl->addPt($pb2, $y2);
                $pl->addPt($pe2, $y2);
            }
            $gd->filledPolygon($pl, $co_pe);
        }
        
        my $fo = "$do/$id.png";
        open(FH, ">$fo") or die "cannot open $fo for write\n";
        print FH $gd->png();
        close FH;
    }
}



