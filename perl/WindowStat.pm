package WindowStat;
use strict;
use InitPath;
use Common;
use Seq;
use Readfile;
use Path::Class;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/write_windows
    window_stats_gc window_stats_gene window_stats_cov window_stats_features
    window_stats_loc
    merge_stats/;
@EXPORT_OK = qw//;
sub window_stats_gc {
    my ($fi, $fo, $refDb) = rearrange(['fi', 'fo', 'refdb'], @_);
    my $t = readTable(-in=>$fi, -header=>1);
    my $f_genome = "$DIR_Genome/$refDb/41_genome.fa";
    my $h = readSeq($f_genome, "fasta", 2);
    my $chrLenH = { map {$_=>getSeqLen($_, $refDb)} uniq($t->col("chr")) };
    open(FH, ">$fo");
    print FH join("\t", qw/chr beg end bp_N bp_GC/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($chr, $s, $e) = map {$t->elm($i, $_)} qw/chr beg end/;
        $e = $chrLenH->{$chr} if $e > $chrLenH->{$chr};
        my $seqStr = $h->{$chr}->subseq($s, $e);
        my ($cnt, $cntN) = (0, 0);
        while($seqStr =~ /([GCN])/ig) {
            if($1 eq "N") {
                $cntN ++;
            } else {
                $cnt ++;
            }
        }
        print FH join("\t", $chr, $s, $e, $cntN, $cnt)."\n";
    }
}
sub window_stats_cov {
    my ($fi, $fo, $p) = rearrange(['fi', 'fo', 'p'], @_);
    my $t = readTable(-in=>$fi, -header=>1);
    my $chrLenH = { map {$_=>getSeqLen($_, $p->{refdb})} uniq($t->col("chr")) };
    open(FH, ">$fo");
    print FH join("\t", qw/chr beg end bp_covered bp_covered_uniq/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($chr, $s, $e) = map {$t->elm($i, $_)} qw/chr beg end/;
        $e = $chrLenH->{$chr} if $e > $chrLenH->{$chr};
        my $loc = Bio::Location::Simple->new(-seq_id=>$chr, -start=>$s, -end=>$e, -strand=>1);
        my ($cnt1) = getPctCov(-loc=>$loc, -acc=>$p->{accs}, -covopt=>1, -p=>$p->{cutoff});
        my ($cnt2) = getPctCov(-loc=>$loc, -acc=>$p->{accs}, -covopt=>2, -p=>$p->{cutoff});
        print FH join("\t", $chr, $s, $e, $cnt1, $cnt2)."\n";
    }
}
sub window_stats_gene {
    my ($fi, $fo, $refDb, $opt) = rearrange(['fi', 'fo', 'refdb', 'opt'], @_);
    my $t = readTable(-in=>$fi, -header=>1);
    my $ld = Localdb->new(-db=>$refDb);
    my $chrLenH = { map {$_=>seqLen($_, $refDb)} uniq($t->col("chr")) };
    open(FH, ">$fo");
    print FH join("\t", qw/chr beg end bp_cds bp_intron bp_utr5 bp_utr3 bp_intergenic id source type note/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($chr, $s, $e) = map {$t->elm($i, $_)} qw/chr beg end/;
        $e = $chrLenH->{$chr} if $e > $chrLenH->{$chr};
        my $loc = Bio::Location::Simple->new(-seq_id=>$chr, -start=>$s, -end=>$e, -strand=>1);
        my $pos = [[$s, $e]];
        my $len = $e - $s + 1;

        my @fes = $ld->getFeatures(-loc=>$loc, -types=>'mRNA');
        my @poss_c = map {[$_->start, $_->end]} $ld->getFeatures(-loc=>$loc, -types=>'CDS');
        my @poss_5 = map {[$_->start, $_->end]} $ld->getFeatures(-loc=>$loc, -types=>'five_prime_UTR');
        my @poss_3 = map {[$_->start, $_->end]} $ld->getFeatures(-loc=>$loc, -types=>'three_prime_UTR');
        my @poss_m = map {[$_->start, $_->end]} @fes;
        my @poss = (@poss_c, @poss_5, @poss_3, @poss_m); 
        my @types = ((1)x@poss_c, (2)x@poss_5, (3)x@poss_3, (4)x@poss_m);

        my ($posCs, $pos5s, $pos3s, $posIs) = ([], [], [], []);
        if(@poss_m) {
            my $ref = tiling(\@poss, \@types, 1);
            for (@$ref) {
                my ($s, $e, $idx) = @$_;
                if($types[$idx] == 1) {
                    push @$posCs, [$s, $e];
                } elsif($types[$idx] == 2) {
                    push @$pos5s, [$s, $e];
                } elsif($types[$idx] == 3) {
                    push @$pos3s, [$s, $e];
                } else {
                    die "type error: unknown type $types[$idx]\n" unless $types[$idx] == 4;
                    push @$posIs, [$s, $e];
                }
            }
        }

        my $cntC = @$posCs ? [posOvlp($pos, $posCs)]->[1] : 0;
        my $cntI = @$posIs ? [posOvlp($pos, $posIs)]->[1] : 0;
        my $cnt5 = @$pos5s ? [posOvlp($pos, $pos5s)]->[1] : 0;
        my $cnt3 = @$pos3s ? [posOvlp($pos, $pos3s)]->[1] : 0;
        my $cntT = @poss_m ? [posCmp($pos, \@poss_m)]->[4] : $len;

        die join(" ", $cntC, $cntI, $cnt5, $cnt3, $cntT, "!=", $len)." at $chr:$s-$e\n" if $cntC + $cntI + $cnt5 + $cnt3 + $cntT != $len;
        my ($id, $source, $type, $note) = ("") x 4;
        if(@fes > 0) {
            my (@ids, @sources, @types, @notes);
            for my $fe (@fes) {
                my @descs = $fe->get_tag_values("Note");
                my $desc = @descs ? $descs[0] : "";
                push @notes, $desc;
                push @ids, $fe->id;
                push @sources, $fe->source;
                push @types, $fe->primary_tag;
            }
            $id = join(" ", @ids);
            $type = join(" ", @types);
            $note = join(" ", @notes);
            $source = join(" ", @sources);
        }
        print FH join("\t", $chr, $s, $e, $cntC, $cntI, $cnt5, $cnt3, $cntT, $id, $source, $type, $note)."\n";
    }
}
sub window_stats_features {
    my ($fi, $fo, $refDb, $types) = rearrange(['fi', 'fo', 'refdb', 'types'], @_);
    $types = [$types] unless ref($types) eq 'ARRAY';
    my $t = readTable(-in=>$fi, -header=>1);
    my $ld = Localdb->new(-db=>$refDb);
    my $chrLenH = { map {$_=>getSeqLen($_, $refDb)} uniq($t->col("chr")) };
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/chr beg end/, @$types)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($chr, $s, $e) = map {$t->elm($i, $_)} qw/chr beg end/;
        $e = $chrLenH->{$chr} if $e > $chrLenH->{$chr};
        my $loc = Bio::Location::Simple->new(-seq_id=>$chr, -start=>$s, -end=>$e, -strand=>1);
        my @lens;
        for my $type (@$types) {
            my @fes = $ld->getFeatures(-loc=>$loc, -types=>$type);
            my @locs = map {[$_->start, $_->end]} @fes;
            my $len = 0;
            if(@locs) {
                my $locs2 = posMerge(\@locs);
                for (sort {$a->[0] <=> $b->[0]} @$locs2) {
                    my ($gS, $gE) = @$_;
                    $gS = $s if $gS < $s;
                    $gE = $e if $gE > $e;
                    $len += $gE - $gS + 1;
                    die "illegal position: $gS - $gE \n" if $len < 0;
                }
            }
            push @lens, $len;
        }
        print $fh join("\t", $chr, $s, $e, @lens)."\n";
    }
}
sub window_stats_loc {
    my ($fi, $fw, $fo, $f_seq) = rearrange([qw/fi fw fo f_seq/], @_);
    my $ti = readTable(-in=>$fi, -header=>1);
    my $tw = readTable(-in=>$fw, -header=>1);

    my @chrs = uniq($tw->col("chr"));
    my $chrLenH = { map {$_=>seqLen($_, $f_seq)} @chrs };

    my @types = uniq($ti->col("cat"));
    print "calculating density for these types:\n\t".join(" ", @types)."\n";
    my $idx_cat = first_index {$_ eq "cat"} $ti->header;
    my $idx_chr = first_index {$_ eq "chr"} $ti->header;
    my $locH;
    for my $type (@types) {
        for my $chr (@chrs) {
            my $t2 = $ti->match_pattern("\$_->[$idx_cat] eq '$type' && \$_->[$idx_chr] eq '$chr'");
            my $loc = [];
            if($t2->nofRow > 0) {
                $loc = [ map {[$t2->elm($_, "beg"), $t2->elm($_, "end")]} (0..$t2->nofRow-1) ];
            }
            $locH->{$type}->{$chr} = $loc;
        }
    }
    
    my @types_cnt = map {"cnt_".$_} @types;
    open(FH, ">$fo");
    print FH join("\t", qw/chr beg end/, @types, @types_cnt)."\n";
    for my $i (0..$tw->nofRow-1) {
        my ($chr, $b, $e) = map {$tw->elm($i, $_)} qw/chr beg end/;
        $e = $chrLenH->{$chr} if $e > $chrLenH->{$chr};
        my $locW = [[$b, $e]];
        my @lens;
        my @cnts;
        for my $type (@types) {
            my $locT = $locH->{$type}->{$chr};
            my ($locO, $len) = posOvlp($locW, $locT);
            push @lens, $len;
            push @cnts, scalar(@$locO);
        }
        print FH join("\t", $chr, $b, $e, @lens, @cnts)."\n";
    }
}


sub merge_stats {
    my ($fis, $fo, $opt) = rearrange(['ins', 'out', 'opt'], @_);
    $opt ||= 1;
    my $t;
    for my $i (0..@$fis-1) {
        my $fi = $fis->[$i];
        my $t2 = readTable(-in=>$fi, -header=>1);
        if($i == 0) {
            $t = $t2;
        } else {
            for my $col ($t2->header) {
                next if $col =~ /(chr)|(start)|(beg)|(end)/;
                $t->addCol($t2->colRef($col), $col);
            }
        }
    }
    open(FH, ">$fo");
    print FH $t->tsv(1);
    close FH;
    if($opt == 2) {
        for my $chr (uniq($t->col("chr"))) {
            my $t4 = $t->match_pattern("\$_->[0] eq '$chr'");
            my $fo2 = $fo;
            $fo2 =~ s/(\.\w+)$/\_$chr$1/;
            open(FH, ">$fo2");
            print FH $t4->tsv(1);
            close FH;
        }
    }
}




1;
__END__
