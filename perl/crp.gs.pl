#!/usr/bin/perl -w
use strict;
use lib ($ENV{"SCRIPT_HOME_PERL"});
use InitPath;
use Common;
use Seq;
use Gtb;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $ver = "20130131";
my $org;
$org = "Athaliana";
#$org = "Mtruncatula_3.5";
#$org = "Osativa";
my $dir = "$DIR_Misc2/crp.gs/$ver/$org";
my $f_ref = "$DIR_Misc4/genome/$org/01_refseq.fa";

my $f01 = "$dir/01.gtb";
my $f02 = "$dir/02_curation.tbl";
my $f_ext = "$DIR_Misc4/spada.crp.$org.simple/31_model_evaluation/61_final.gtb";
my $f04 = "$dir/04_updated.gtb";
my $f11 = "$dir/11_new.gtb";
#update_gs($f01, $f_ext, $f02, $f04, $f11);
#gtb2Gff($f11, "$dir/11_new.gff");
sub update_gs {
    my ($fg, $fge, $fa, $fo, $fon) = @_;
    my $tg = readTable(-in=>$fg, -header=>1);
    my $ta = readTable(-in=>$fa, -header=>1);

    my $tge = readTable(-in=>$fge, -header=>1);
    my $he = { map {$tge->elm($_, "id") => $_} (0..$tge->nofRow-1) };

    $ta = $ta->match_pattern("\$_->[-1] eq 1");
    my @idxs = map {$he->{$_}} $ta->col("id");
    my $tn = $tge->subTable(\@idxs, [$tge->header]);
    printf "%s: %d new models added\n", $org, $tn->nofRow;

    for my $i (0..$tn->nofRow-1) {
        $tg->addRow( $tn->rowRef($i) );
    }
    open(FH, ">$fo") or die "cannot open $fo for writing\n";
    print FH $tg->tsv(1);
    close FH;
    open(FH, ">$fon") or die "cannot open $fon for writing\n";
    print FH $tn->tsv(1);
    close FH;
}

#$f04 = $f01;
my $f05 = "$dir/05_renamed.gtb";
#crp_rename($f04, $f_ref, $f05);
sub crp_rename {
    my ($fi, $f_ref, $fo) = @_;
    my $ti = readTable(-in=>$fi, -header=>1);
    $ti->sort("cat3", 1, 0, "chr", 1, 0, "beg", 0, 0);
    
    my @chrs = uniq($ti->col("chr"));
    my @len_digits = map {getDigits(seqLen($_, $f_ref) / 1000000)} @chrs;
    my $chr_digits = getDigits(scalar(grep /\d+/, @chrs));
    my $h = { map {$chrs[$_] => $len_digits[$_]} 0..$#chrs };

    for my $i (0..$ti->nofRow-1) {
        my ($chr, $beg, $fam) = map {$ti->elm($i, $_)} qw/chr beg cat3/;
        my $begStr = sprintf "%0".$h->{$chr}."d", $beg/1000000;
        my $chrStr = $chr;
        $chrStr =~ s/chr//i;
        $chrStr = sprintf "%0".$chr_digits."d", $chrStr if $chrStr =~ /^\d+$/;

        my $id = sprintf "%s_chr%s_%sM", lc($fam), $chrStr, $begStr;
        $ti->setElm($i, "id", $id);
    }
    
    my $ref = group($ti->colRef("id"));
    my $hd = { map { $_ => getDigits($ref->{$_}->[1]) } keys(%$ref) };
    my $hc;
    for my $i (0..$ti->nofRow-1) {
        my $id = $ti->elm($i, "id");
        $hc->{$id} ||= 0;
        my $cnt = ++$hc->{$id};
        
        $id = sprintf "%s_%0".$hd->{$id}."d", $id, $cnt;
        $ti->setElm($i, "id", $id);
        my $pa = "$id.p";
        $ti->setElm($i, "parent", $pa);
    }
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO $ti->tsv(1);
    close FHO;
}
#gtb2Gff($f05, "$dir/05.gff");

my $d10 = "$dir/10_sigp";
my $f10_01 = "$d10/01.gtb";
my $f10_06 = "$d10/06.fa";
my $f10_11 = "$d10/11_sigp.tbl";
#crp_stat_sigp($f10_01, $f10_06, $f10_11);

sub crp_stat_sigp {
    my ($f_gtb, $f_seq, $fo) = @_;
    my $tg = readTable(-in=>$f_gtb, -header=>1);
    my $h_seq = readSeq($f_seq, 2);
  
    open(FH, ">$fo");
    print FH join("\t", qw/id parent cat2 cat3 note sigp sp_pos/)."\n";
    for my $i (0..$tg->nofRow-1) {
        my ($id, $pa, $chr, $beg, $end, $strand, $locE, $locI, $locC, $loc5, $loc3, $phase, $source, $conf, $cat1, $cat2, $cat3, $note) = $tg->row($i);
        my $seq_pro = $h_seq->{$pa}->seq;
        die "no sequence for $id\n" unless $seq_pro;
        $seq_pro =~ s/\*$//;
        my $seq = Bio::Seq->new(-id=>$id, -seq=>$seq_pro);
        my ($sp, $sp_pos) = run_signalp($seq, 3);
        print FH join("\t", $id, $pa, $cat2, $cat3, $note, $sp, $sp_pos)."\n";
    }
    close FH;
}
sub run_signalp {
    my ($seq, $version) = @_;
    my $fTmp = "sigp.fa";
    writeSeq($seq, $fTmp);
    my ($prob1, $prob2, $site) = ("") x 3;
    open(OUT, "signalp -t euk $fTmp |") or die "failed: $!\n";
    my ($pos, $score, $sp) = ("", "", 0);
    if($version == 3) {
        my ($score1, $score2);
        while( <OUT> ) {
            chomp;
            if(/Signal peptide probability: ([\d\.]+)/i) {
                $score1 = $1;
            } elsif(/Signal anchor probability: ([\d\.]+)/i) {
                $score2 = $1;
            } elsif(/Max cleavage site probability: [\d\.]+ between pos\. (\d+)/i) {
                $pos = $1+1;
            }
        }
        $score = max($score1, $score2);
        $sp = ($score >= 0.450) ? 1 : 0;
        $pos = "" if $sp == 0;
    } elsif($version == 4) {
        while( <OUT> ) {
            chomp;
            next unless /^tmp/;
            my ($id, $Cmax, $posC, $Ymax, $posY, $Smax, $posS, $Smean, $D, $tag, $Dmaxcut, $network) = split " ";
            ($pos, $score, $sp) = ($posY, $D, 1) if $tag eq "Y";
        }
    } else {
        die "unsupported version: $version\n";
    }
    system("rm $fTmp");
    return ($sp, $pos);
}



