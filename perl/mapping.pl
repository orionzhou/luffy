#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use File::Path qw/make_path remove_tree/;
use InitPath;
use Mapping;
use Data::Dumper;
use Gff;
use Gtb;
my @pKs = qw/fedb refdb best gap pctidty pctcov/;
my $pVs = {
    11 => [qw/mt_35_gi   mt_35 1 3000 0.9/],
    12 => [qw/mt_35_gea  mt_35 1 3000 0.9/],
    13 => [qw/mt_35_defl mt_35 1 3000 0.9/],
    14 => [qw/mt_35_nbs  mt_35 1 3000 0.9/],
    15 => [qw/mt_35_ncr  mt_35 1 3000 0.9/],
    16 => [qw/mt_35_crp  mt_35 1 3000 0.9/],
    21 => [qw/mt_gi      mt_30 1 3000 0.9/],
    22 => [qw/mt_gea     mt_30 1 3000 0.9/],
    23 => [qw/mt_defl    mt_30 1 3000 0.9/],
    31 => [qw/mt_30      mt_35 1 1000 0.9 0.1/],
    32 => [qw/mt_35_map  mt_30 1 3000 0.8/],
    41 => [qw/mt_35_cgh  mt_35 1 10   0.9/],
    42 => [qw/mt_35_affy mt_35 1 1    0.9/],
    43 => [qw/hm340_affy hm340 1 1    0.9/],
    51 => [qw/gm_snp     Gmax  1 10   0.8/],
    62 => [qw/crp_at     Athaliana  1 3000 0.9/],
    64 => [qw/crp_os     Osativa  1 3000 0.9/],
};
my $opt = 62;
my $p = { map { $pKs[$_] => $pVs->{$opt}->[$_] } 0..$#pKs };
my $dir = "$DIR_misc2/mapping/$opt\_".$p->{fedb};
make_path($dir) unless -d $dir;

my $f01 = "$dir/01_seq.fa";

my $f_tgt = "$DIR_db/blat/".$p->{refdb}.".2bit";
pipe_blat(-qry=>$f01, -tgt=>$f_tgt, -dir=>$dir, -p=>$p);

#pipe_gmap(-fseq=>$f01, -dir=>$dir, -p=>$p);

#pipe_manual($dir);

#$p->{ld} = Localdb->new(-db=>$p->{fedb}, -opt=>$p->{opt});
#$p->{ld}->loadGff(-empty=>1, -files=>$f15);

__END__
my $f21 = file($dir, "21_sum.txt");
#changeFormat(-in=>$f05, -out=>$f06, -db1=>$p->{fedb}, -db2=>$p->{refdb});
#printBacSum($f06);
#mappingStat($f06, "mt_30_id", "mt_35_loc");
my $f22 = file($dir, "22_refined.txt");
#refineMapping($f21, $f22);
#printBacSum($f22);
my $f23 = file($dir, "23_sum.txt");
#addGene(-in=>$f07, -out=>$f08, -db=>$p->{refdb});
#geneSum($f08);
my $f23x = file($dir, "23x_sum.txt");
#convLoc(-in=>$f23, -out=>$f23x, -refdb=>$p->{refdb}, -opt=>3);
sub changeFormat {
    my ($fi, $fo, $db1, $db2) = rearrange(['in', 'out', 'db1', 'db2'], @_);
    my $t = readTable(-in=>$fi, -header=>1);
    $t->header([qw/mt_30_id mt_35_chr mt_35_loc pct_cov pct_idty mt_30_loc/]);
    my $tmp = $t->delCol("mt_30_loc");
    $t->addCol($tmp, "mt_30_loc", 1);
    $tmp = $t->delCol("mt_35_chr");
    $t->addCol($tmp, "mt_35_chr");
    $tmp = $t->delCol("mt_35_loc");
    $t->addCol($tmp, "mt_35_loc");
    $t->addCol([('') x $t->nofRow], "mt_30_chr", 1);
    $t->addCol([('') x $t->nofRow], "bac_support");
    my ($ld1, $ld2) = map {Localdb->new(-db=>$_)} ($db1, $db2);
    for my $i (0..$t->nofRow-1) {
        my @ps = split(" ", $t->elm($i, "mt_30_loc"));
        my ($mt_30_chr, $mt_30_loc) = split(":", $ps[1]);
        $t->setElm($i, "mt_30_chr", $mt_30_chr);
        $t->setElm($i, "mt_30_loc", $mt_30_loc);
        next unless $t->elm($i, "mt_35_chr");
        my ($mt_35_chr, $mt_35_loc) = map {$t->elm($i, $_)} qw/mt_35_chr mt_35_loc/;
        my $loc1G = locStr2Obj($mt_30_loc, $mt_30_chr);
        my $loc1L = $ld1->chr2Bac($loc1G);
        my $bac1 = getMainSeqId($loc1L);
        die "invalid mt_30_loc [$mt_30_loc]\n".Dumper($loc1L) if $bac1 =~ /^gap/i;
        my $loc2G = locStr2Obj($mt_35_loc, $mt_35_chr);
        my $loc2L = $ld2->chr2Bac($loc2G);
        my $bac2 = getMainSeqId($loc2L);
        die "invalid mt_35_loc [$mt_35_chr $mt_35_loc]\n".Dumper($loc2L) if $bac2 =~ /^gap/i ;
        my $bacSupport = $bac1 eq $bac2 ? 1 : "$bac1 $bac2";
        $t->setElm($i, "bac_support", $bacSupport);
    }
    my $fh = new IO::File $fo, "w";
    print $fh $t->tsv(1);
}
sub printBacSum {
    my ($fi) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $idx35 = first_index {$_ eq "mt_35_loc"} $t->header;
    my $cnt1 = $t->nofRow;
    $t = $t->match_pattern("\$_->[$idx35] ne ''"); 
    my $cnt2 = $t->nofRow;
    my ($cnt3, $cnt4, $cnt5) = (0) x 3;
    my ($cnt6, $cnt7, $cnt8) = (0) x 3;
    print "a total of $cnt1 lines, ".($cnt1-$cnt2)." are unmapped\n";
    my $stat = { map {$_=>$t->elm($_, "bac_support") eq 1 ? 1 : 0} (0..$t->nofRow-1) };
    my $ref = group($t->colRef("mt_30_id"));
    print scalar(keys %$ref)." ids are mapped to ".$t->nofRow." locations\n";
    my @ids1 = grep {$ref->{$_}->[1] == 1} keys %$ref;
    my @idxs1 = map {$ref->{$_}->[0]} @ids1;
    my @idxs1b = grep {$stat->{$_} eq 1} @idxs1;
    printf "\t%5d ids are uniquely mapped\n", scalar(@ids1);
    printf "\t\t%5d are consistent\n", scalar(@idxs1b);
    my @ids2 = grep {$ref->{$_}->[1] > 1} keys %$ref;
    printf "\t%5d ids are non-uniquely mapped to %5d locations\n", scalar(@ids2), sum(map {$ref->{$_}->[1]} @ids2);
    for my $id (@ids2) {
        my ($idx, $n) = @{$ref->{$id}};
        my $stat2 = { map {$_ => $stat->{$_}} ($idx..$idx+$n-1) };
        my $tag = sum(values %$stat2);
        if($tag == 0) {
            $cnt3 ++;
            $cnt6 += $n;
        } elsif($tag == 1) {
            $cnt4 ++;
            $cnt7 += $n;
        } else {
            die "fatal error (H)\n" unless $tag > 1;
            $cnt5 ++;
            $cnt8 += $n;
        }
    }
    printf "\t\t%5d ids have  0 consistent mappings[%5d]\n", $cnt3, $cnt6;
    printf "\t\t%5d ids have  1 consistent mappings[%5d]\n", $cnt4, $cnt7;
    printf "\t\t%5d ids have >1 consistent mappings[%5d]\n", $cnt5, $cnt8;
}
sub refineMapping {
    my ($fi, $fo) = rearrange(['in', 'out'], @_);
    my $t = readTable(-in=>$fi, -header=>1);
    my $stat = { map {$_=>$t->elm($_, "bac_support") eq 1 ? 1 : 0} (0..$t->nofRow-1) };
    my $ref = group($t->colRef("mt_30_id"));
    my @idxsRm;
    my @ids2 = grep {$ref->{$_}->[1] > 1} keys %$ref;
    for my $id (@ids2) {
        my ($idx, $n) = @{$ref->{$id}};
        my $stat2 = { map {$_ => $stat->{$_}} ($idx..$idx+$n-1) };
        for my $i ($idx..$idx+$n-1) {
            die "not id: $id\n" unless $t->elm($i, "mt_30_id") eq $id;
        }
        if( sum(values %$stat2) >= 1 ) {
            for my $idxRm (grep {$stat->{$_} == 0} keys %$stat2) {
                push @idxsRm, $idxRm;
                delete $stat2->{$idxRm};
            }
        }
    }
    print "\t".@idxsRm." inconsistent mappings removed\n";
    $t->delRows(\@idxsRm);
    my $fh = new IO::File $fo, "w";
    print $fh $t->tsv(1);
}
sub addGene {
    my ($fi, $fo, $db) = rearrange(['in', 'out', 'db'], @_);
    my $t = readTable(-in=>$fi, -header=>1);
    my $ld = Localdb->new(-db=>$db);
    my @buf;
    for my $i (0..$t->nofRow-1) {
        unless($t->elm($i, "mt_35_loc")) {push @buf, ""; next;};
        my $loc = locStr2Obj(map {$t->elm($i, $_)} qw/mt_35_loc mt_35_chr/);
        my $seqLen = $loc->length;
        my @fes = $ld->getFeatures(-loc=>$loc, -types=>['gene', 'transposable_element_gene']);
        my @tmp;
        for my $fe (@fes) {
            my $feLen = $fe->length;
            push @tmp, $fe->id."[$feLen/$seqLen]";
        }
        push @buf, join(" ", @tmp);
    }
    $t->addCol(\@buf, "mt_35_id");
    my $fh = new IO::File $fo, "w";
    print $fh $t->tsv(1);
}
sub geneSum {
    my ($fIn) = @_;
    my $t = readTable(-in=>$fIn, -header=>1);
    my $idx1 = first_index {$_ eq "mt_35_loc"} $t->header;
    my $t2 = $t->match_pattern("\$_->[$idx1] ne ''");
    my $cntL = $t2->nofRow;
    my @genes = $t2->col("mt_35_id");
    my $cntE = grep {!$_} @genes;
    my $cntM = grep /\s/, @genes;
    my $cntU = grep /^\S+$/, @genes;
    printf "\t%5d locations:\n", $cntL;
    printf "\t\t%5d have  0 gene(s)\n", $cntE;
    printf "\t\t%5d have  1 gene(s)\n", $cntU;
    printf "\t\t%5d have >1 gene(s)\n", $cntM;
}



