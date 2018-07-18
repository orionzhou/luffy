package Writefile;
use strict; use Init; use Common; use Localdb; use Readfile; use Seq; use Parser;
use Path::Class; use Run; use IO::File; use Data::Table; use LWP::Simple; use URI::Escape;
use Bio::AlignIO; use Bio::Seq; use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/translate mergeFile mergeChr getFeDetail
    checkDnaSeq writeSeq6Rf writeCatR writeCatDefl writeSsp
    extractIdByRegex writeCovById convLoc
    writeVntSum writeInfoById fileRegex sortFileByCol processExp/;
@EXPORT_OK = qw//;
sub translate {
    #translate(-in=>$fCDS, -out1=>$fProtein, -out2=>$fInfo);
    my ($fIn, $fOut1, $fOut2) = rearrange(['in', 'out1', 'out2'], @_);
    my $fOutH = new IO::File $fOut2, "w";
    my $seqInH = Bio::SeqIO->new(-file => "<".$fIn, -format => "fasta");
    my $seqOutH = Bio::SeqIO->new(-file => ">".$fOut1, -format => 'fasta');
    my @lenAry;
    my ($cnt, $lenThresh) = (0, 50);
    while(my $seq = $seqInH->next_seq()) {
        my ($noteStart, $noteStop) = ("", "");
        my $lenExp = int($seq->length/3);
        if($seq->seq !~ /^atg/i) {
            my $lenBefore = $seq->length;
            my $seqStr = $seq->seq;
            $seqStr =~ s/^.*?(atg.*)/$1/i;
            my $lenAfter = length($seqStr);
            if($lenBefore == $lenAfter) {
                $noteStart = "No start codon";
                $seq = Bio::Seq->new(-id=>$seq->id, -desc=>$seq->description, -seq=>"");
            } else {
                $noteStart = "start-> ".($lenBefore-$lenAfter)." bp";
                $seq = Bio::Seq->new(-id=>$seq->id, -desc=>$seq->description, -seq=>$seqStr);
            }
        }
        my $prot = $seq->translate(-complete => 1);
        my $protStr = $prot->seq;
        if($prot->seq =~ /^(.*)\*/) {
            $protStr = $1;
            my $lenBefore = $prot->length;
            my $lenAfter = length($protStr);
            $noteStop = "stop-> $lenAfter/$lenBefore aa";
            $prot = Bio::Seq->new(-id=>$prot->id, -desc=>$prot->description, -seq=>$protStr);
        }
        my ($lenObs) = ($prot->length);
        push (@lenAry, $lenObs);
        if($lenObs>$lenThresh) {
            $seqOutH->write_seq($prot);
            $cnt ++;
        }
        print $fOutH join("\t", $seq->id, $lenExp, $lenObs, $noteStart, $noteStop, $prot->description, "\n");
    }
    print $cnt." protein sequences passed filtering\n";
}
sub getFeDetail {
    my ($fo, $ids, $refDb) = rearrange(['out', 'ids', 'refdb'], @_);
    my $fh = new IO::File $fo, "w";
    my $ld = Localdb->new(-db=>$refDb);
    for my $id (@$ids) {
        my $fe = $ld->getFeatureByName($id);
        print $fh join("\t", $id, $fe->seq_id, $fe->start, $fe->end, $fe->source, join(" ", $fe->get_tag_values("Note")))."\n";
    }
}
sub writeSeq6Rf {
#my $f00 = file($DIR_Genome, "mt_35", '41_genome.fa');
#my $f01 = file($DIR_Misc2, "hmmsearch", "02_crp_hmmsearchX", '01_seq.fa');
#writeSeq6Rf(-ins=>[$f00], -out=>$f01);
    my ($fis, $fo) = rearrange(['ins', 'out'], @_);
    my ($size, $step) = (300_000) x 2;
    my $seqOH = Bio::SeqIO->new(-file=>">$fo", -format=>'fasta');
    for my $fi (@$fis) {
        my $seqIH = Bio::SeqIO->new(-file=>$fi, -format=>'fasta');
        while(my $seq= $seqIH->next_seq) {
            my $seqLen = $seq->length;
            my $it = getWindowIter(1, $seqLen, $size, $step);
            while(my $ref = $it->()) {
                my ($s, $e) = @$ref;
                my $seqStrF = $seq->subseq($s, $e);
                my $seqStrR = Bio::Seq->new(-seq=>$seqStrF)->revcom->seq;
                for my $strand ('+', '-') {
                    my $seqStr = $strand eq "-" ? $seqStrR : $seqStrF;
                    for my $rf (0..2) {
                        my $id2 = join("_", $seq->id, $s, $e, $strand, $rf);
                        my $seqDna  = Bio::Seq->new(-id=>$id2, -seq=>$seqStr, -alphabet=>'dna');
                        my $seqPro = $seqDna->translate(-frame=>$rf);
                        $seqOH->write_seq($seqPro);
                        print join("\t", $id2)."\n";
                    }
                }
            }
        }
    }
}
sub writeCatR {
    my ($fIn, $fOut) = rearrange(['in', 'out'], @_);
    $fIn ||= file($DIR_In, "gene_id.txt");
    $fOut ||= file($DIR_In, "gene_id_r.txt");
    my $fOutH = new IO::File $fOut, "w";
    print $fOutH join("\t", "family", "gene")."\n";
    my $catHash = getId(-in=>$fIn, -format=>1);
    for my $cat (sort keys %$catHash) {
        my $idHash = $catHash->{$cat};
        for my $id (sort keys %$idHash) {
            print $fOutH join("\t", $cat, $id)."\n";
        }
    }
}
sub writeCatDefl {
    my ($fIn, $fOut) = rearrange(['in', 'out'], @_);
    $fIn ||= file($DIR_In, 'defl_cat.txt');
    $fOut ||= file($DIR_In, 'id_defl.txt');
    my $fInH = new IO::File $fIn, "r";
    my $fOutH = new IO::File $fOut, "w";
    my $catHash = getId(-in=>$fIn, -format=>1);
    my ($biodb) = getDBI('mt_defl'); 
    my @feAry = getFeatures(-type=>['mRNA'], -db=>$biodb);
    my $fe2Cat = {};
    for my $fe (@feAry) {
        next if $fe->seq_id !~ /MtChr/i;
        my $note = join(" |", $fe->get_tag_values("Note"));
        while(my ($cat, $idHash) = each(%$catHash)) {
            for my $id2 (keys %$idHash) {
                if($note =~ /$id2/i) {
                    $fe2Cat->{$fe->id} = $cat;
                    last;
                }
            }
        }
        $fe2Cat->{$fe->id} = "un-categorized" unless exists $fe2Cat->{$fe->id};
    }
    my $tmp = {};
    for my $id (sort {$fe2Cat->{$a} cmp $fe2Cat->{$b} || $a cmp $b} keys %$fe2Cat) {
        if(!exists($tmp->{$fe2Cat->{$id}})) {
            print $fOutH ">".$fe2Cat->{$id}."\n";
            $tmp->{$fe2Cat->{$id}} = 1;
        }
        print $fOutH join("\t", $id)."\n";
    }
}
sub writeGmSnp {
    my ($fIn, $param) = @_;
    my ($ref1, $ref2) = readTable(-in=>$fIn, -head=>$param->{head}, -left=>$param->{left});
    my (@pos, $htAcc);
    for my $ref (@$ref2) {
        push @pos, $ref->[-1];
    }
    for my $acc (sort(keys(%$ref1))) {
        my $ref = $ref1->{$acc};
        my @snpAry = @$ref[$param->{left}-1..scalar(@$ref)-1];
        die "not ".scalar(@pos)." SNPs at {$acc}[".scalar(@snpAry)."]\n" unless @snpAry == @pos;
        my @htAry;
        for my $ele (@snpAry) {
            my @pair;
            $ele = uc($ele);
            $ele =~ s/\?/N/g;
            if(length($ele) == 3) {
                @pair = split("/", $ele);
            } elsif(length($ele) == 2) {
                @pair = split("", $ele);
            } elsif(length($ele) == 1) {
                @pair = ($ele) x 2;
            }
            for my $i (0..$#pair) {
                my $snp = $pair[$i];
                if($snp !~ /^[ATCGN]$/) {
                    print "unknown allele[$snp] at $acc converted to N\n";
                    $snp =~ s/[^ATCGN]/N/;
                }
                $htAry[$i] .= $snp;
            }
        }
        $htAcc->{$acc} = [ @htAry ];
    }
    print "\tSNP [".scalar(@pos)."]; Acc [".scalar(keys(%$htAcc))."]\n";
    $fIn =~ s/\.\w+$//;
    my $myOut = VntOut->new(-pos=>[@pos], -seq=>$htAcc);
    $myOut->write(-format=>"haploview", -out=>$fIn);
}
sub extractIdByRegex {
#my ($cat, $regex) = qw/NCR late nodulin|cysteine\-rich/;
#($cat, $regex) = qw/NB-ARC NB\-ARC|NB\-LRR|NBS\-LRR/;
#my $f01 = file($DIR_Out, "id_$family.txt");
#extractIdByRegex(-regex=>$regex, -cat=>$cat, -out=>$f01, -db=>);
    my ($regex, $fOut, $cat, $db) = rearrange(['regex', 'out', 'cat', "db"], @_);
    my $fOutH = new IO::File $fOut, "w";
    my $ld = Localdb->new(-db=>$db);
    my @fe = $ld->getFeatures(-type=>'mRNA');
    my $rst = {};
    for my $fe (@fe) {
        next if $fe->seq_id =~ /chr0/i;
        my $note = join(" |", $fe->get_tag_values("Note"));
        $rst->{$fe->id} = [$fe->source, $note] if $note =~ /$regex/i;
    }
    print $fOutH ">$cat\n";
    for my $feId (sort(keys(%$rst))) {
        print $fOutH join("\t", $feId, @{$rst->{$feId}})."\n";
    }
}
sub writeCovById {
    my ($idAry, $accAry, $fOut, $feDb, $cutOffCov) = 
        rearrange(['ids', 'acc', 'out', 'fedb', 'cutoffcov'], @_);
    my $fOutH = new IO::File $fOut, "w";
    $cutOffCov ||= 2;
    my $cnt = 0;
    print $fOutH join("\t", qw/id length acc avg_Cov avg_UniqCov pct_Cov pct_UniqCov/)."\n";
    my $ld = Localdb->new(-db=>$feDb);
    for my $id (sort @$idAry) {
        my $fe = $ld->getFeatureByName($id);
        my $len = $fe->end - $fe->start + 1;
        my $loc = $fe->seq_id.":".$fe->start."..".$fe->end;
        my $covH1 = getCov(-acc=>$accAry, -loc=>$loc, -covopt=>1);
        my $covH2 = getCov(-acc=>$accAry, -loc=>$loc, -covopt=>2);
        for my $acc (sort @$accAry) {
            my @covAry1 = @{$covH1->{$acc}};
            my @covAry2 = @{$covH2->{$acc}};
            my $avgCov1 = sprintf("%.02f", sum(@covAry1) / $len);
            my $avgCov2 = sprintf("%.02f", sum(@covAry2) / $len);
            my $cntCov1 = sum( map {$_ >= $cutOffCov ? 1 : 0} @covAry1 );
            my $cntCov2 = sum( map {$_ >= $cutOffCov ? 1 : 0} @covAry2 );
            my $pctCov1 = sprintf("%.03f", $cntCov1/$len);
            my $pctCov2 = sprintf("%.03f", $cntCov2/$len);
            print $fOutH join("\t", $id, $len, $acc, $avgCov1, $avgCov2, $pctCov1, $pctCov2)."\n";
        }
        print $cnt." done\n" if ++$cnt % 100 == 0;
    }
}
sub writeVntSum {
    my ($fId, $fVnt, $fCov, $fOut, $accAry) =
        rearrange(['fid', 'fvnt', 'fcov', 'out', 'acc'], @_);
    my $catH = getId(-in=>$fId, -format=>1);
    my $idHash = getId(-in=>$fId, -format=>2);
    my @idAry = keys(%$idHash);
    my $vntSum = {};
    for my $id (@idAry) {
        for my $acc (@$accAry) {
            $vntSum->{$id}->{$acc} = [(0) x 6] unless exists $vntSum->{$id}->{$acc};
        }
    }
    my $fVntH  = new IO::File $fVnt, "r";
    my $fCovH  = new IO::File $fCov, "r";
    my $fOutH  = new IO::File $fOut, "w";
    print $fOutH join("\t", qw/cat id acc length avgCov avgUniqCov pctCov pctUniqCov vntCnt/)."\n";
    while( <$fVntH> ) {
        chomp;
        next if /^\#/;
        my @ps = split($_);
        my ($gene, $acc) = @ps[0, 2];
        die "$gene - $acc not exists\n" unless exists $vntSum->{$gene}->{$acc};
        my $ref = $vntSum->{$gene}->{$acc};
        @$ref[0..4] = @ps[1, 3..6];
    }
    while( <$fCovH> ) {
        chomp;
        next if /^\#/;
        my @ps = split($_);
        my ($gene, $acc) = @ps[0,4];
        die "$gene - $acc not exists\n" unless exists $vntSum->{$gene}->{$acc};
        $vntSum->{$gene}->{$acc}->[5] ++;
    }
    for my $cat (sort(keys(%$catH))) {
        for my $id (sort(keys(%{$catH->{$cat}}))) {
            die "$id not in vntSum\n" unless exists $vntSum->{$id};
            my $statAcc = $vntSum->{$id};
            for my $acc (sort(keys(%$statAcc))) {
                print $fOutH join("\t", $cat, $id, $acc, @{$statAcc->{$acc}})."\n";
            }
        }
    }
}
sub writeInfoById {
    my ($fIn, $fOut, $col, $feDb, $opt, $tag) = 
        rearrange(['in', 'out', 'col', 'fedb', 'opt', 'tag'], @_);
    my @ids = getId(-in=>$fIn, -col=>$col); 
    my $fOutH = new IO::File $fOut, "w";
    my $ld = Localdb->new(-db=>$feDb);
    print $fOutH join("\t", qw/id source desc/)."\n" if $opt == 1;
    for my $id (@ids) {
        my $fe = $ld->getFeatureByName($id);
        die "error in $id\n" unless $fe;
        if($opt == 1) {
            print $fOutH join("\t", $id, $fe->source, join("|", $fe->get_tag_values("Note")))."\n";
        } elsif($opt == 2) {
            print $fOutH join("\t", $fe->seq_id, $fe->start, $tag)."\n";
        } else {
            die "unsupported opt[$opt]\n";
        }
    }
}
sub fileRegex {
    my ($fIn, $fOut, $find, $replace) = @_;
    die "$fIn is not there\n" unless -s $fIn;
    my $fInH = new IO::File $fIn, "r";
    my $fOutH = new IO::File $fOut, "w";
    while( <$fInH> ) {
        my $line = $_;
        $line =~ s/\Q$find/$replace/g;
        print $fOutH $line;
    }
}
sub sortFileByCol {
#sortFileByCol(-in=>$fIn, -out=>$fOut, -cols=>['3n','4n'])
    my ($fIn, $fOut, $colAry) = rearrange(['in', 'out', 'cols'], @_);
    die "$fIn is not there\n" unless -s $fIn;
    $colAry = ref($colAry) ? $colAry : [$colAry];
    my $colStr = join(" ", map {"-k".$_} @$colAry);
    my $cmd = qq/sort $colStr -S 90% -o $fOut $fIn/;
    open(JJ, $cmd." |") || die "Failed: $! in \n$cmd\n";
    while( <JJ> ) {
        print "\t".$_."\n";
    }
}
sub processExp {
    my ($fIn, $fOut) = rearrange(['in', 'out'], @_);
    my $fInH = new IO::File $fIn, "r";
    my $fOutH = new IO::File $fOut, "w";
    my $cnt = 0;
    while( <$fInH> ) {
        $_ =~ s/\r?\n$//;
        if( /^\#/ ) {
            print $fOutH $_."\n";
            next;
        }
        my @ps = split("\t");
        if( $ps[0] =~ /contig\_(\d+)\_fgenesh\_gene(\d+)/ ) {
            $ps[0] = "ctg".$1."_".int($2).".1";
        }
        for my $i (1..$#ps) {
            $ps[$i] = sprintf("%.03f", $ps[$i]);
        }
        my $str = join("\t", @ps);
        print $fOutH "$str\n";
        last if ++$cnt < 1;
    }
}
sub convLoc {
    my ($fi, $fo, $db) = rearrange(['in', 'out', 'refdb'], @_);
    my $ld = Localdb->new(-db=>$db);
    open(IN, "<$fi");
    open(OU, ">$fo");
    while(<IN>) {
        chomp;
        if(/^(\#)|(chr)/) {
            print OU $_."\n";
            next;
        }
        my @ps = split "\t";
        my ($c, $str, $s, $e) = @ps[0,3,4,6];
        if($c =~ /^chr[0UT]$/) {
            my $loc1 = Bio::Location::Simple->new(-seq_id=>$c, -strand=>$str, -start=>$s, -end=>$e);
            my $loc2 = $ld->chr2Bac($loc1);
            die "hits on >1 bacs:\n".Dumper($loc1).Dumper($loc2) unless $loc2->is_single_sequence();
            @ps[0,3,4,6] = map {$loc2->$_} qw/chr strand beg end/;
        } else {
            die "unknown seq_id: $c\n" unless $c =~ /^chr([1-8])$/;
        }
        print OU join("\t", @ps)."\n";
    }
}




1;
__END__
