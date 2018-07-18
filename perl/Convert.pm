package Convert;
use strict; use Init; use Common; use Path::Class; use IO::File; use File::Copy;
use Bio::AlignIO; use Bio::Seq; use Bio::SeqIO; use Data::Dumper; use Seq;
use VntWrite; use Readfile; use Parser;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/psl2Mtb bls2Mtb mtb2Htb pgp2Tbl
    tbl2Tlf
    convertTree convertAln convertSeq fastphase2haploview ms2MigrateN/;
@EXPORT_OK = qw//;
sub convertAln {
#convertAln(-indir=>$dirIn, -outdir=>$dirOut, -informat=>'fasta', -outformat=>'phylip');
    my ($dirIn, $dirOut, $inFormat, $outFormat) = rearrange(
        ['indir', 'outdir', 'informat', 'outformat'], @_);
    my $formatSuffix = {"phylip"=>"phy", "fasta"=>"fa", "clustalw"=>"aln"};
    die("unsupported input  format $inFormat\n")  unless exists($formatSuffix->{$inFormat});
    die("unsupported output format $outFormat\n") unless exists $formatSuffix->{$outFormat};
    system('mkdir', $dirOut) unless -d $dirOut;
    opendir(IN, $dirIn) || print "Can't open $dirIn";
    my @dirContent = readdir(IN);
    my $cntFile = 0;
    for my $child (@dirContent) {
        if($child !~ /^\./) {
            if(-d dir($dirIn, $child)) {
                print "$child\t";
                convertAln(dir($dirIn, $child), dir($dirOut, $child), $inFormat, $outFormat);
            } elsif(-s file($dirIn, $child)) {
                my $fIn  = file($dirIn, $child);
                my ($inSuffix, $outSuffix) = ($formatSuffix->{$inFormat}, $formatSuffix->{$outFormat});
                my $outChild = $child;
                $outChild =~ s/\.$inSuffix$/\.$outSuffix/;
                my $fOut = file($dirOut, $outChild);
                my $in   = Bio::AlignIO->new(-file => "<".$fIn,  -format => 'fasta');
                my $out  = Bio::AlignIO->new(-file => ">".$fOut, -format => 'phylip' );
                while ( my $aln = $in->next_aln() ) {
                    $out->write_aln($aln);
                }
                last if ++$cntFile > 10000000;
            }
        }
    }
    print "$cntFile files converted\n";
}
sub convertTree {
    my ($fi, $fo, $inF, $outF) = rearrange(['in', 'out', 'inf', 'outf'], @_);
    my $treeI = Bio::TreeIO->new(-format=>$inF, -file=>$fi);
    my $tree = $treeI->next_tree;
    die("no tree in there\n") unless $tree;
    my $treeO = Bio::TreeIO->new(-format=>$outF, -file=>'>'.$fo);
    $treeO->bootstrap_style('traditional') if $outF eq 'newick';
    for my $node ( $tree->get_leaf_nodes ) {
        my $id = $node->id;
        my @anno;
    }
    $treeO->write_tree($tree);
}
sub convertSeq {
    my ($fi, $fo, $inF, $outF) = rearrange(['in', 'out', 'inf', 'outf'], @_);
    my $seqs = readSeq($fi, $inF);
    writeSeq(-seqs=>$seqs, -out=>$fo, -format=>$outF);
}
sub checkHetero {
    my ($ref) = @_;
    my @accs = sort keys %$ref;
    my @seqs = map {$ref->{$_}} @accs;
    my @cnts;
    for my $acc (@accs) {
        my ($ht1, $ht2) = @{$ref->{$acc}};
        my $cnt = 0;
        for my $i (0..@$ht1-1) {
            my ($nt1, $nt2) = ($ht1->[$i], $ht2->[$i]);
            $cnt ++ if $nt1 ne $nt2;
            $ht2->[$i] = $ht1->[$i];
        }
        push @cnts, $cnt;
    }
    return (\@accs, \@seqs, \@cnts);
}
sub fastphase2haploview {
    my ($inDir1, $inDir2, $fWin, $outDir) = rearrange(["indir1", "indir2", 'fwin', "outdir"], @_);
    die "inDir1 is not there\n" unless -d $inDir1;
    die "inDir2 is not there\n" unless -d $inDir2;
    system("mkdir -p $outDir") unless -d $outDir;
    my $fo = file($outDir, "..", "fastphase_hetero.txt");
    my $foh = new IO::File $fo, "w";
    my $accs;
    my $t = readTable(-in=>$fWin, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($round, $cntSnp) = map {$t->elm($i, $_)} qw/round cntSnp/;
        next unless $cntSnp > 1;
        last if $i < 0;
        my $fi1 = file($inDir1, "$round\_hapguess_indiv.out");
        my $fi2 = file($inDir2, "$round\_fastphase.txt");
        die "$fi1 is not there\n" unless -s $fi1;
        die "$fi2 is not there\n" unless -s $fi2;
        my $ref = parse_fastphase($fi1);
        my ($accs2, $seqs, $stats) = checkHetero($ref);
        unless($accs) {
            $accs = $accs2;
            print $foh join("\t", qw/round cntSnp/, @$accs)."\n";
        } else {
            my ($aryS, $ary1, $ary2) = aryCmp($accs, $accs2);
            die "acc not uniform:\n".Dumper($accs).Dumper($accs2) if @$ary1 || @$ary2;
        }
        my $fh = new IO::File $fi2;
        my $line = readline($fh);
        $line = readline($fh);
        $line = readline($fh);
        my @poss = split(" ", $line);
        shift @poss;
        my $vnt = VntOut->new(-pos=>\@poss, -acc=>$accs, -seq=>$seqs);
        $vnt->write(-pre=>file($outDir, $round), -formats=>"haploview");
        print $foh join("\t", $round, $cntSnp, @$stats)."\n";
    }
}
sub ms2MigrateN {
    my ($fi, $fo) = @_;
    my $fhi = new IO::File $fi, "r";
    my $fho = new IO::File $fo, "w";
    my @stats = split(" ", readline($fhi));
    my $idx = first_index {$_ eq "-I"} @stats;
    die "no subpop in $fi\n" unless $idx >= 0;
    my $nRep = $stats[2];
    my $nPop = $stats[$idx+1];
    my @pSizes = @stats[$idx+2..$idx+1+$nPop];
    my $tSize = sum(@pSizes);
    print "$nRep reps:\n";
    print "\t$nPop pops: ".join("+", @pSizes)."\n";
    my @data;
    my @tmp;
    my $tag = 0;
    while(<$fhi>) {
        chomp;
        if(!$_) {
            $tag = 0;
            push @data, [@tmp] if @tmp;
            @tmp = ();
        } elsif(/^\/\//) {
            $tag = 1;
            next;
        }
        if($tag == 1) {
            push @tmp, $_;
        }
    }
    push @data, [@tmp] if @tmp;
    die "not $nRep reps: ".@data."\n" unless $nRep == @data;
    print $fho join("\t", $nPop, $nRep, "[ms]")."\n";
    my @lens = map {[split(" ", $_->[0])]->[1]} @data;
    print $fho join("\t", @lens)."\n";
    my $idxP = 2;
    for my $i (0..$nPop-1) {
        my $pSize = $pSizes[$i];
        print $fho join("\t", $pSize, "Pop$i")."\n";
        for my $j (0..$nRep-1) {
            my @rows = @{$data[$j]}[$idxP..$idxP+$pSize-1];
            for my $k (0..$#rows) {
                my $str = $rows[$k];
                $str =~ tr/01/AT/;
                printf $fho "$i\_%02d\t%s\n", $k, $str;
            }
        }
        $idxP += $pSize;
    }
}

sub tbl2Tlf {
#my $dir = dir($DIR_Repo, "mt_35/31_phylogeny");
#my $fi = file($dir, "mt_label.tbl");
#my $fo = file($dir, "mt_label.tlf");
#tbl2Tlf($fi, $fo);
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $fh = new IO::File $fo, "w";
    my @colnames = $t->header;
    my $idx_id = 0;
    for my $i (0..$t->nofRow-1) {
        my @ps = $t->row($i);
        my $id = $ps[$idx_id];
        print $fh join(" ", $id, map {$colnames[$_]." {".$ps[$_]."}"} (0..$#colnames))."\n";
    }
    close $fh;
}




1;
__END__
