#!/usr/bin/perl
use strict;
use lib $ENV{"SCRIPT_HOME_PERL"};
use Data::Dumper;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use Time::HiRes qw/gettimeofday tv_interval/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use LWP::Simple;
use Web;
require LWP::Parallel::UserAgent;
use HTTP::Request; 

use InitPath;
use Common;
use Seq;

my $org = "Athaliana";
#$org = "Osativa";
my $dir = "$DIR_Misc2/crp.ssp/$org";
my $f05 = "$dir/05.tbl";

#transform_tair9("$dir/01.tbl", $f05, "$dir/../tair9_updates.tbl");

my $f01 = "$dir/01_ids.txt";
#get_oryspssp($f01, $f05);
#get_oryspssp_p($f01, $f05);

my $f_ref = "$DIR_Misc4/genome/$org/01_refseq.fa";
my $f06 = "$dir/06_seq.tbl";
#prepare_dna_seq_os($f05, $f06, $f_ref);
#prepare_dna_seq_at($f05, $f06, $f_ref);
my $f07 = "$dir/07_rel.tbl";
#run_genewise($f06, $f07);
my $f08 = "$dir/08.gtb";
#convert2Gtb($f07, $f08, $f_ref);

sub map_pos {
    my ($pos, $ary, $opt) = @_;
    my $idx = first_index {$_->[0] >= $pos} @$ary;
    
    my $posN;
    if( $idx == -1 ) {
        my ($beg1, $beg2) = @{$ary->[-1]};
        die "unknown pos [$pos] > $beg1\n" if $pos < $beg1;
        $posN = $pos - $beg1 + $beg2;
    } else {
        my ($end1, $end2) = @{$ary->[$idx]};
        if($pos == $end1) {
            $posN = $end2;
        } else {
            my ($beg1, $beg2) = @{$ary->[$idx-1]};
            if( $end2 - $beg2 == 1 ) {
                $posN = $opt eq "from" ? $end2 : $beg2;
            } elsif( ($end1-$beg1) == ($end2-$beg2) ) {
                $posN = $pos - $beg1 + $beg2;
            } else {
                die "cannot map $pos [$beg1-$end1] [$beg2-$end2]\n";
            }
        }
    }
    return $posN;
}
sub transform_tair9 {
    my ($fi, $fo, $fc) = @_;

    my $h;
    my $tc = readTable(-in=>$fc, -header=>1);
    for my $i (1..$tc->nofRow) {
        my ($chr, $pos1, $pos2) = $tc->row($i-1);
        $h->{$chr} ||= [];
        push @{$h->{$chr}}, [$pos1, $pos2];
    }

    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO join("\t", qw/id chr srd beg end seq_pro/)."\n";
    
    my $t = readTable(-in=>$fi, -header=>1);
    for my $i (1..$t->nofRow) {
        my ($id, $chr, $srd, $beg, $end, $seq_sig, $seq_pro_pro) = $t->row($i-1);
        
        $chr = "Chr$chr";
        ($beg, $end) = ($end, $beg) if $beg > $end;
        $srd = $srd eq "bottom" ? "-" : "+";
        $seq_pro_pro =~ s/\s+$//;
        my $seq_pro = $seq_sig.$seq_pro_pro;

        $beg = map_pos($beg, $h->{$chr}, "from");
        $end = map_pos($end, $h->{$chr}, "to");

        print FHO join("\t", $id, $chr, $srd, $beg, $end, $seq_pro)."\n";
    }
    close FHO;
}
sub get_oryspssp {
    my ($fi, $fo) = @_;

    my @ids;
    open(FHI, "<$fi") or die "cannot open $fi for reading\n";
    while(<FHI>) {
        chomp;
        push @ids, $_;
    }
    close FHI;
    print "\t".@ids." IDs read\n";

    my @keys = (qw/source chromosome strand start stop dnaseq preprotein/, 
        "signalpeptide", "target organelle", "nr annotation",
        "domain annotation", "gene in the locus", "left neighbor gene", "right neighbor gene"
    );
    
    my $t0 = [gettimeofday];
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO join("\t", "id", @keys)."\n";
    for my $i (0..$#ids) {
        my $id = $ids[$i];
        my $url_base = "http://www.genoportal.org/PSSP/showGeneDetail.do";
        my $url = "$url_base?name=$id";
        my $str = get($url);
        
        my $h = { map {$_=>""} @keys };
        my $te = HTML::TableExtract->new(depth=>1, count=>1);
        $te->parse($str);
        for my $ts ($te->tables) {
            for my $row ($ts->rows) {
                $h->{$row->[0]} = $row->[1] if exists $h->{$row->[0]};
            }
        }

        my @values = map {$h->{$_}} @keys;
        print FHO join("\t", $id, @values)."\n"; 
        if( ($i+1) % 1000 == 0 ) {
            printf "   %5d: %.01f min\n", $i+1, tv_interval($t0, [gettimeofday]) / 60;
        }
    }
    close FHO;
}
sub get_oryspssp_p {
    my ($fi, $fo) = @_;

    my @ids;
    open(FHI, "<$fi") or die "cannot open $fi for reading\n";
    while(<FHI>) {
        chomp;
        push @ids, $_;
    }
    close FHI;
    print "\t".@ids." IDs read\n";

    my @keys = (qw/source chromosome strand start stop dnaseq preprotein/, 
        "signalpeptide", "target organelle", "nr annotation",
        "domain annotation", "gene in the locus", "left neighbor gene", "right neighbor gene"
    );

    my $url_base = "http://www.genoportal.org/PSSP/showGeneDetail.do";
    
    my $total = @ids;
    my $per_cnt = 1000;
    my $cnt = ceil($total/$per_cnt);

    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO join("\t", "id", @keys)."\n";
    
    for my $i (1..$cnt) {
        my $idx_beg = $per_cnt * ($i-1);
        my $idx_end = min($total-1, $per_cnt*$i-1);
        
        my $reqs = [ map {HTTP::Request->new('GET', "$url_base?name=$_")} @ids[$idx_beg..$idx_end] ],

        my $pua = LWP::Parallel::UserAgent->new();
        $pua->in_order  (1);
        $pua->duplicates(1);
        $pua->timeout   (200);
        $pua->redirect  (1);

        my $t0 = [gettimeofday];
        foreach my $req (@$reqs) {
            if ( my $res = $pua->register ($req) ) { 
                print STDERR $res->error_as_HTML; 
            }  
        }
        my $entries = $pua->wait();

        foreach (keys %$entries) {
            my $res = $entries->{$_}->response;
            $res->request->url =~ /name=(\w+)/;
            my $id = $1;
            my $str = $res->content;
            
            my $h = { map {$_=>""} @keys };
            my $te = HTML::TableExtract->new(depth=>1, count=>1);
            $te->parse($str);
            for my $ts ($te->tables) {
                for my $row ($ts->rows) {
                    $h->{$row->[0]} = $row->[1] if exists $h->{$row->[0]};
                }
            }

            my @values = map {$h->{$_}} @keys;
            print FHO join("\t", $id, @values)."\n"; 
        }
        printf "   %5d - %5d: %.01f min\n", $idx_beg+1, $idx_end+1, tv_interval($t0, [gettimeofday]) / 60;
        $t0 = [gettimeofday];
    }
    close FHO;
}
sub prepare_dna_seq_os {
    my ($fi, $fo, $f_ref) = @_;
    
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO join("\t", qw/id src chr begr endr srd loc seq_pro seq/)."\n";
    
    my $hl;
    my $relax = 12;
    my $ti = readTable(-in=>$fi, -header=>1);
    for my $i (1..$ti->nofRow) {
        my ($id, $src, $chr, $srd, $beg, $end, $seq_dna_s, $seq_pro) = $ti->row($i-1);
        die join("\t", $ti->row($i))."\n" if $chr eq "";
        $chr = "Chr$chr";
        $hl->{$chr} ||= seqLen($chr, $f_ref);
        
        if($end > $hl->{$chr}) {
            print "$id: $chr $beg-$end ignored\n";
            next;
        }
        my $begr = max(1, $beg - $relax);
        my $endr = min($end + $relax, $hl->{$chr});
        my $seq_dna = seqRet([[$begr, $endr]], $chr, $srd, $f_ref);
        print FHO join("\t", $id, $src, $chr, $begr, $endr, $srd, '', $seq_pro, $seq_dna)."\n";
        printf "  %5d / %5d done\r", $i, $ti->nofRow;
    }
    print "\n";
    close FHO;
}
sub prepare_dna_seq_at {
    my ($fi, $fo, $f_ref) = @_;
    
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO join("\t", qw/id src chr begr endr srd loc seq_pro seq/)."\n";
    
    my $hl;
    my $relax = 12;
    my $t = readTable(-in=>$fi, -header=>1);
    for my $i (1..$t->nofRow) {
        my ($id, $chr, $srd, $beg, $end, $seq_pro) = $t->row($i-1);
        
        $hl->{$chr} ||= seqLen($chr, $f_ref);
        if($end > $hl->{$chr}) {
            print "$id: $chr $beg-$end ignored\n";
            next;
        }

        my $begr = max(1, $beg - $relax);
        my $endr = min($end + $relax, $hl->{$chr});
        my $seq_dna = seqRet([[$begr, $endr]], $chr, $srd, $f_ref);
        print FHO join("\t", $id, '', $chr, $begr, $endr, $srd, '', $seq_pro, $seq_dna)."\n";
        printf "  %5d / %5d done\r", $i, $t->nofRow;
    }
    print "\n";
    close FHO;
}
sub diff_pro {
    my ($str1, $str2) = @_;
    my ($len1, $len2) = (length($str1), length($str2));
    my $diff = 0;
    for my $i (0..(min($len1, $len2)-1)) {
        $diff ++ if substr($str1, $i, 1) ne substr($str2, $i, 1);
    }
    $diff += abs($len1 - $len2);
    return $diff;
}
sub sum_genewise1 {
    my ($fi) = @_;
    open(FHT, "<$fi") || die "cannot open $fi for reading\n";
    my $loc = [];
    while(<FHT>) {
        chomp;
        my @ps = split "\t";
        if(@ps == 9 && $ps[1] eq "GeneWise" && $ps[2] eq "cds") {
            push @$loc, [$ps[3], $ps[4]];
        } 
    }
    close FHT;
    return $loc;
}
sub run_genewise {
    my ($fi, $fo) = @_;
    
    $ENV{"GeneWise"} = "/home/msi/zhoup/spada_soft/wise2.2.0";
    my $f_bin = $ENV{"GeneWise"}."/bin/genewise";
    $ENV{"WISECONFIGDIR"} = $ENV{"GeneWise"}."/wisecfg";
    die "$f_bin not there\n" unless -s $f_bin;
    
    $ENV{"TMP_DIR"} ||= ".";
    my $f_dna = $ENV{"TMP_DIR"}."/genewise_dna.fa";
    my $f_pro = $ENV{"TMP_DIR"}."/genewise_pro.fa";
    my $f_out = $ENV{"TMP_DIR"}."/genewise.txt";
 
    my $ti = readTable(-in=>$fi, -header=>1);
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO join("\t", $ti->header)."\n";
    for my $i (1..$ti->nofRow) {
        my ($id, $src, $chr, $begr, $endr, $srd, $locS, $seq_pro_o, $seq) = $ti->row($i-1);
        next if $i < 0;
        
        writeFile($f_dna, ">seq", $seq);
        writeFile($f_pro, ">pro", $seq_pro_o);
        runCmd("$f_bin -gff $f_pro $f_dna > $f_out", 0);
        my $loc = sum_genewise1($f_out);

        if( @$loc == 0 ) {
#            print join("\t", $id, $chr, $begr, $endr, $srd)."\n";
        } else {
            $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
            my $seq_dna = getSubSeq($seq, $loc);
            my $seq_pro = Bio::Seq->new(-seq=>$seq_dna)->translate->seq;
            
            if($seq_pro !~ /\*$/) {
                my $seq_end = getSubSeq($seq, [[$loc->[-1]->[1]+1, $loc->[-1]->[1]+3]]);
                if($seq_end =~ /^(TAA)|(TAG)|(TGA)$/ ) {
                    $loc->[-1]->[1] += 3;
                    $seq_pro .= "*";
                }
            }
            
            $seq_pro_o .= "*" if $seq_pro eq "$seq_pro_o\*";
            if(diff_pro($seq_pro, $seq_pro_o) < 3) {
                $ti->setElm($i-1, "loc", locAry2Str($loc));
            } else {
                print join("\t", $id, $src, $chr, $begr, $endr, $srd, $seq_pro, $seq_pro_o)."\n";
            }
        }
        print FHO join("\t", $ti->row($i-1))."\n";
        printf "  %5d / %5d done\r", $i, $ti->nofRow;
    }
    print "\n";
    system("rm $f_out $f_dna $f_pro");
    close FHO;
}
sub convert2Gtb {
    my ($fi, $fo, $f_ref) = @_;
    
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO join("\t", qw/id parent chr beg end strand locE locI locC loc5 loc3 phase source conf cat1 cat2 cat3 note seq/)."\n";

    my $ti = readTable(-in=>$fi, -header=>1);
    my $t = $ti->match_pattern("\$_->[6] ne ''");
    printf " %d out of %d passed\n", $t->nofRow, $ti->nofRow;
    for my $i (1..$t->nofRow) {
        my ($id, $src, $chr, $begr, $endr, $srd, $locS, $seq_pro_o, $seq) = $t->row($i-1);
        
        my $locL = locStr2Ary($locS);
        $locL = [sort {$a->[0] <=> $b->[0]} @$locL];
        my $loc = $srd eq "-" ? [ reverse map {[$endr-$_->[1]+1, $endr-$_->[0]+1]} @$locL ]
            : [ map {[$_->[0]+$begr-1, $_->[1]+$begr-1]} @$locL ];
        my $seq_dna = seqRet($loc, $chr, $srd, $f_ref);
        my $seq_pro = Bio::Seq->new(-seq=>$seq_dna)->translate->seq;
        if(diff_pro($seq_pro_o, $seq_pro) >= 3) {
            print join("\t", $id, $src, $chr, $srd, $locS)."\n";
            print "$seq_pro\n";
            print "$seq_pro_o\n";
            die;
        }

        $loc = [sort {$a->[0] <=> $b->[0]} @$loc];
        my ($beg, $end) = ($loc->[0]->[0], $loc->[-1]->[1]);
        my $locM = [[$beg, $end]];
        my ($locI) = posDiff($locM, $loc);
        my $phaseS = join(",", getPhase($loc, $srd));
        print FHO join("\t", $id, "$id.p", $chr, $beg, $end, $srd, 
            locAry2Str($loc), locAry2Str($locI), locAry2Str($loc), '', '', $phaseS, 
            $src, '', 'gene', 'mRNA', '', '', $seq_pro )."\n";
        printf "  %5d / %5d done\r", $i, $t->nofRow;
    }
    print "\n";
    close FHO;
}





