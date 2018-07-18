package Localdb;
use strict; use Init; use Common; use Path::Class; use LWP::UserAgent;
use Bio::DB::SeqFeature::Store; use Bio::AlignIO; use Bio::Seq; 
use Data::Dumper; use Seq; use DBI;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Clone qw/clone/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/getBioDb getDbh 
    getExp getWindows partWindows
    assignWindows getLongestOrf extendLoc getGeneStructure/;
@EXPORT_OK = qw/getFeatures getFeatureByName getNote 
    loadGff seqRetById get_transcript_from_gene_id 
    getRegionStat getLoc retLeadSeq bac2Chr chr2Bac/;
sub new {
    my $self = shift;
    my ($dbname, $opt) = rearrange(['db', 'opt'], @_);
    $dbname ||= "mt_35";
    my ($dbi, $dbh, $biodb, $type, $subtype);
    my @db_mthap_list = qw/mt_30 mt_35 gm_101/;
    my $idxDb = first_index {$_ eq $dbname} @db_mthap_list;
    if($idxDb != -1) {
        $dbname = 'mt_hapmap' if $dbname eq "mt_30";
        $dbi = join(":", "dbi", "mysql", $dbname, $ENV{"MYSQLH_1"});
        $dbh = DBI->connect($dbi, $ENV{"MYSQLU_1"}, $ENV{"MYSQLP_1"}, {"RaiseError"=>1});
        $biodb = Bio::DB::SeqFeature::Store->new(
            -dsn=>$dbi, -write=>1, -user=>$ENV{'MYSQLU_1'}, -pass=>$ENV{'MYSQLP_1'})
            or die "Cannot connect to [$dbname]";
        $opt ||= "mRNA";
    } else {
        $dbi = join(":", "dbi", "mysql", $dbname, $ENV{"MYSQLH_2"});
        $dbh = DBI->connect($dbi, $ENV{"MYSQLU_2"}, $ENV{"MYSQLP_2"}, {"RaiseError"=>1});
        $biodb = Bio::DB::SeqFeature::Store->new(
            -dsn=>$dbi, -write=>1, -user=>$ENV{'MYSQLU_2'}, -pass=>$ENV{'MYSQLP_2'})
            or die "Cannot connect to [$dbname]";
        $opt ||= "match";
    }
    if($opt eq "gene") {
        ($type, $subtype) = ('gene', 'mRNA');
    } elsif($opt eq "mRNA") {
        ($type, $subtype) = ('mRNA', [qw/exon CDS five_prime_UTR three_prime_UTR/]);
    } elsif($opt eq "match") {
        ($type, $subtype) = ('match', 'match_part');
    } else {
        die("unknown option[$opt]\n");
    }
    return bless {dbh=>$dbh, biodb=>$biodb, type=>$type, subtype=>$subtype}, 
        ref($self) || $self;
}
sub getDbh {
    my ($db) = @_;
    unless( $db ) {
        $db = "mt_35";
        print "refDb not specified - using $db as default\n";
    }
    my $dbname = $db;
    if($db eq "mt_30") {
        $dbname = "mt_hapmap";
    }
    my $dbi = join(":", "dbi", "mysql", $dbname, $ENV{"MYSQLH_1"});
    my $dbh = DBI->connect($dbi, $ENV{"MYSQLU_1"}, $ENV{"MYSQLP_1"}, {"RaiseError"=>1});
    return $dbh;
}
sub getBioDb {
    my ($db) = @_;
    unless( $db ) {
        $db = "mt_35";
        print "refDb not specified - using $db as default\n";
    }
    my $dbname = $db;
    if($db eq "mt_30") {
        $dbname = "mt_hapmap";
    } elsif($db eq "mt_35") {
    } elsif($db eq "gm_101") {
    }
    my $dbi = join(":", "dbi", "mysql", $dbname, $ENV{"MYSQLH_1"});
    my $dbh = DBI->connect($dbi, $ENV{"MYSQLU_1"}, $ENV{"MYSQLP_1"}, {"RaiseError"=>1});
    my $biodb = Bio::DB::SeqFeature::Store->new(-dsn=>$dbi, -write=>1, -user=>$ENV{'MYSQLU_1'}, -pass=>$ENV{'MYSQLP_1'}) or die "Cannot connect to db";
    die "unknown refDb: $db\n" unless $biodb;
    return $biodb;
}
sub getFeatures {
    my $self = shift;
    my ($locObj, $types) = rearrange(['loc', 'types'], @_);
    my $biodb = $self->{biodb};
    my @fes;
    if($locObj){ 
        my $seqid = $locObj->seq_id;
        die "not uniform seqid\n".Dumper($locObj) unless $seqid;
        my @locs = $locObj->each_Location();
        for (@locs) {
            push @fes, $biodb->features(-seqid=>$seqid, -start=>$_->start, -end=>$_->end, -types=>$types);
        }
    } else {
        @fes = $biodb->features(-types=>$types);
    }
    return @fes;
}
sub getFeatureByName {
    my $self = shift;
    my $biodb = $self->{biodb};
    my ($id) = @_;
    my @fes = $biodb->get_features_by_name($id);
    my @note;
    if(@fes == 0) {
        return undef;
    } elsif(@fes > 1) {
        return $fes[0];
    } else { 
        return $fes[0];
    }
}
sub getNote {
    my $self = shift;
    my ($id) = @_;
    my $fe = $self->getFeatureByName($id);
    my @note;
    if($fe->primary_tag eq "gene") {
        my @mRNA = $fe->get_SeqFeatures("mRNA");
        for my $oneMRNA (@mRNA) {
            push (@note, $oneMRNA->get_tag_values("Note"));
        }
    } elsif($fe->primary_tag =~ /(prime_utr|cds|exon|intron)/i) {
        my $parent = $fe->get_tag_values("Parent");
        die "no parent for ".$fe->id."\n" unless $parent;
        my $mRNA = $self->getFeatureByName($parent);
        push (@note, $mRNA->get_tag_values("Note"));
    } else {
        push (@note, $fe->get_tag_values("Note"));
    }
    return join(" |", @note);
}
sub loadGff {
    my $self = shift;
    my ($empty, $fAry) = rearrange(["empty", "files"], @_);
    $empty ||= 0;
    $fAry = [$fAry] unless ref($fAry) eq 'ARRAY';
    for my $fIn (@$fAry) {
        die "$fIn is not there\n" unless -s $fIn;
    }
    my $biodb = $self->{biodb};
    if( $empty == 1) {
        print " !!! database emptyed |||\n";
        $biodb->init_database('erase');
    }
    my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(-store=>$biodb, -verbose=>1)
        or die "Couldn't create GFF3 loader";
    my @fHs;
    for my $fIn (@$fAry) {
        print $fIn."\n";
        my $fInH = new IO::File $fIn, "r";
        push @fHs, $fInH;
    }
    $loader->load(@fHs);
}
sub partWindows {
    my ($locAry, $winSize, $winStep) = rearrange(['loc', 'winsize', 'winstep'], @_);
    $winSize ||= 0;
    $winStep ||= 0;
    $winStep = ($winStep>0 && $winStep<=$winSize) ? $winStep : $winSize;
    my $rst = [];
    $locAry = posMerge($locAry);
    if($winSize > 0) {
        for my $locPair (sort {$a->[0]<=>$b->[0]} (@$locAry)) {
            my ($frameBeg, $frameEnd) = (@$locPair);
            my $round = ceil( ($frameEnd - $frameBeg + 1) / $winStep );
            for my $i (0..($round-1)) {
                my ($beg, $end) = ($frameBeg + $i*$winStep, $frameBeg + $i*$winStep + $winSize-1);
                $end = $end<$frameEnd ? $end : $frameEnd;
                push @$rst, [$beg, $end];
                last if $end == $frameEnd;
            }
        }
    } else {
        $rst = $locAry;
    }
    return $rst;
}
sub getWindows {
    my ($chr, $winSize, $winStep, $db, $opt) = rearrange(["chr", "winsize", "winstep", "db", "opt"], @_);
    $opt ||= 1;
    $chr = "chr$chr" if $chr =~ /^\d$/;
    my $biodb = getBioDb($db);
    my $locAry = [];
    if($opt == 1) {
        my @bacs = $biodb->features(-seqid=>$chr, -types=>['bac']);
        push @$locAry, map { [$_->start, $_->end] } sort {$a->start <=> $b->start} @bacs;
    } elsif($opt == 2) {
        my $chrLen = getSeqLen($chr, $db);
        push @$locAry , [1, $chrLen];
    } elsif($opt == 3) {
        my @bacs = $biodb->features(-seqid=>$chr, -types=>['bac']);
        for (sort {$a->start <=> $b->start} @bacs) {
            next unless $_->score == 3;
            push @$locAry, [$_->start, $_->end];
        }
    } else {
        die "unsupported option: $opt\n";
    }
    return partWindows(-loc=>$locAry, -winsize=>$winSize, -winstep=>$winStep);
}
sub assignWindows {
    my ($param, $pos, $chr, $refDb) = rearrange(['param', 'pos', 'chr', 'refdb'], @_);
    my ($wsize, $wstep, $ssize, $sstep) = map {$param->{$_}} qw/size step snpsize snpstep/;
    my $win1 = getWindows(-chr=>$chr, -db=>$refDb, -winsize=>$wsize, -winstep=>$wstep, -opt=>2);
    my (@win, @idx);
    for my $i (0..@$win1-1) {
        my ($wS, $wE) = @{$win1->[$i]};
        my @idx2 = sort {$a<=>$b} grep {$pos->[$_] >= $wS && $pos->[$_] <= $wE} (0..@$pos-1);
        $win[$i] = [$wS, $wE, 0];
        $idx[$i] = [];
        if(@idx2) {
            my $win2 = partWindows(-loc=>[[$idx2[0], $idx2[-1]]], -winsize=>$ssize, -winstep=>$sstep);
            my @idxS = map {[$_->[0]..$_->[1]]} @$win2;
            $idx[$i] = [ @idxS ];
            $win[$i]->[2] = @idxS;
        }
    }
    return (\@win, \@idx);
}
sub get_transcript_from_gene_id {
    my $self = shift;
    my ($id) = @_;
    my $fe = $self->getFeatureByName($id);
    my $feM = undef;
    if(!$fe) {
    } elsif ($fe->primary_tag eq "gene" || $fe->primary_tag eq "transposable_element_gene") {
        my @feAry2 = $fe->get_SeqFeatures("mRNA");
        my @feAry3 = map { [$_->get_SeqFeatures("CDS")] } @feAry2;
        my @lenAry = map { sum(map {$_->length} @$_) } @feAry3;
        $feM = $feAry2[first_index {$_ == max(@lenAry)} @lenAry];
    } else {
        die "$id is not a gene/TE_gene\n";
    }
    return $feM;
}
sub getLoc {
    my $self = shift;
    my ($id, $opt, $refDb) = rearrange([qw/id opt refdb/], @_);
    $opt ||= 1;
    my $fe = $self->getFeatureByName($id);
    die "$id is not a transcript\n" if !$fe || $fe->primary_tag ne "mRNA";
    
    my $seqid = $fe->seq_id;
    my $strand = $fe->strand;
    my ($posM, $posC, $posE, $posU5, $posU3, $posI) = getGeneStructure($fe);
    my $pos = {mrna=>$posM, cds=>$posC, exon=>$posE, utr5=>$posU5, utr3=>$posU3, intron=>$posI};
    if($opt == 2) {
        my ($posDos, $posUps) = ([], []);
        my $len = 1000;
        my $cbeg = min( map {$_->[0]} @$posC );
        my $cend = max( map {$_->[1]} @$posC );
      
        my ($lbeg, $lend) = ($cbeg - $len, $cbeg - 1);
        $lbeg = max($lbeg, 1);
        my $loc = Bio::Location::Simple->new(-seq_id=>$seqid, -start=>$lbeg, -end=>$lend);
        my @fes = $self->getFeatures(-loc=>$loc, -types=>'CDS');
        $lbeg = max( map { $_->end } @fes ) + 1 if @fes > 0;

        my ($rbeg, $rend) = ($cend + 1, $cend + $len);
        $rend = min($rend, getSeqLen($seqid, $refDb));
        $loc = Bio::Location::Simple->new(-seq_id=>$seqid, -start=>$rbeg, -end=>$rend);
        @fes = $self->getFeatures(-loc=>$loc, -types=>'CDS');
        $rend = min( map {$_->start} @fes ) - 1 if @fes > 0;

        my ($ubeg, $uend, $dbeg, $dend);
        ($ubeg, $uend) = $strand == -1 ? ($rbeg, $rend) : ($lbeg, $lend);
        ($dbeg, $dend) = $strand == -1 ? ($lbeg, $lend) : ($rbeg, $rend);
        $posUps = [[ $ubeg, $uend ]] if $ubeg <= $uend;
        $posDos = [[ $dbeg, $dend ]] if $dbeg <= $dend;
        $pos->{up1k} = $posUps;
        $pos->{dw1k} = $posDos;
    }
    return ($seqid, $strand, $pos);
}
sub seqRetById {
    my $self = shift;
    my ($id, $type, $refDb) = rearrange(['id', 'type', 'refdb'], @_);
    my ($seqid, $strand, $pos) = $self->getLoc(-id=>$id, -opt=>1, -refdb=>$refDb);
    die "unknown sequence type: $type\n" unless exists $pos->{$type};
    my $loc = locAry2Obj($pos->{$type}, $strand, $seqid);
    my $f_refseq = file($DIR_Genome, $refDb, "41_genome.fa");
    my $seqStr = seqRet($loc, $f_refseq);
    return Bio::Seq->new(-id=>$id, -seq=>$seqStr, -description=>$type);
}
sub retLeadSeq {
    my $self = shift;
    my ($id, $refDb, $opt) = rearrange(['id', 'refdb', 'opt'], @_);
    my @types;
    if($opt == 1) {
        @types = qw/up1000 in500/;
    } elsif($opt == 2) {
        @types = qw/up200 in100/;
    } else {
        die "unsupported option($opt)\n";
    }
    my $seq1 = $self->seqRetById(-id=>$id, -type=>$types[0], -refdb=>$refDb);
    my $seq2 = $self->seqRetById(-id=>$id, -type=>$types[1], -refdb=>$refDb);
    my $seqDesc = "LeadingSeq";
    $seq1 ||= Bio::Seq->new(-seq=>"", -alphabet=>"dna");
    $seq2 ||= Bio::Seq->new(-seq=>"", -alphabet=>"dna");
    my ($len1, $len2) = ($seq1->length, $seq2->length);
    my $start1 = $len1 ? 1 : 0;
    $seqDesc .= "$start1..$len1 ".($len1+1)."..".($len1+$len2);
    my $seq = Bio::Seq->new(-seq=>$seq1->seq.$seq2->seq, -id=>$id, -description=>$seqDesc);
    return ($seq, $len1, $len2);
}
sub chr2Bac {
    my $self = shift;
    my ($loc1) = @_;
    my @fes = $self->getFeatures(-loc=>$loc1, -types=>['bac', 'contig', 'gap']);
    my $loc2 = Bio::Location::Split->new();
    my ($s1, $e1, $str1, $c1) = map {$loc1->$_} qw/start end strand seq_id/;
    if(@fes == 0) {
        die "no bacs found for $loc1\n";
    } else {
        for my $fe (sort {$a->start <=> $b->start} @fes) {
            my ($c3, $s3, $e3, $str3) = map {$fe->$_} qw/id start end strand/;
            $str3 ||= 1;
            my ($c2, $str2) = ($c3, $str1 * $str3);
            my ($s2, $e2);
            if($s3 > $s1) {
                $s2 = 1 if $str3 == 1;
                $e2 = $e3 - $s3 + 1 if $str3 == -1;
            } else {
                $s2 = $s1 - $s3 + 1 if $str3 == 1;
                $e2 = $e3 - $s1 + 1 if $str3 == -1;
            }
            if($e3 < $e1) {
                $e2 = $e3 - $s3 + 1 if $str3 == 1;
                $s2 = 1 if $str3 == -1;
            } else {
                $e2 = $e1 - $s3 + 1 if $str3 == 1;
                $s2 = $e3 - $e1 + 1 if $str3 == -1;
            }
            ($s2, $e2) = ($e2, $s2) if $s2 > $e2;
            my $subLoc = Bio::Location::Simple->new(-seq_id=>$c2, -start=>$s2, -end=>$e2, -strand=>$str2);
            $loc2->add_sub_Location($subLoc);
        }
    }
    if( $loc2->is_single_sequence ) {
        $loc2->seq_id([$loc2->sub_Location()]->[0]->seq_id);
    }
    return $loc2;
}
sub bac2Chr {
    my $self = shift;
    my ($loc1) = @_;
    my $seqid1 = $loc1->seq_id;
    my $fe = $self->getFeatureByName($seqid1);
    unless($fe) {
        $seqid1 =~ s/\.\d+$//;
        $fe = $self->getFeatureByName($seqid1);
    }
    return undef unless $fe;
    my $s2 = l2G($fe->start, $fe->end, $fe->strand, $loc1->start);
    my $e2 = l2G($fe->start, $fe->end, $fe->strand, $loc1->end);
    ($s2, $e2) = ($e2, $s2) if $s2 > $e2;
    my $loc2 = Bio::Location::Simple->new(-seq_id=>$fe->seq_id, -start=>$s2, -end=>$e2, -strand=>$loc1->strand * $fe->strand );
    return $loc2;
}
sub getFirstLine {
    my ($fi) = @_;
    die "$fi is not there\n" unless -s $fi;
    my $fh = new IO::File $fi, "r";
    my $line = readline($fh);
    $line =~ s/\r?\n$//;
    return $line;
}
sub getExp {
    my ($ids, $ds) = @_;
    $ids = ref($ids) eq "ARRAY" ? $ids : [$ids];
    my $dbh = getBioDb("mt_30");
    my $table = "exp_$ds";
    my $qry = "SELECT value FROM $table WHERE id=?";
    my $sth = $dbh->prepare($qry);
    my ($ref, @condS, @condD);
    my $fi = file($DIR_Misc2, $ds, "02_exp.txt");
    my @ps = split("\t", getFirstLine($fi));
    if($ds eq "mtgea") {
        for my $cond (@ps[1..$#ps]) {
            if($cond =~ /^(\S+) \((\S+)\)/) {
                push @condS, $1;
                push @condD, $2;
            } else {
                push @condS, $cond;
                push @condD, "";
            }
        }
    } elsif($ds eq "rnaseq") {
        @condS = @ps[1..$#ps];
    }
    my $cntHit = 0;
    for my $id ( @$ids ) {
        $sth->execute($id) or die "Couldn't execute query: " . $sth->errstr;
        die $sth->rows." hits for $id\n" unless $sth->rows <= 1;
        if( $sth->rows == 0 ) {
            print "\tno exp data for $id\n";
        } else {
            $cntHit ++;
            my $exp = $sth->fetchrow_hashref()->{value};
            $ref->{$id} = [ split(" ", $exp) ];
        }
    }
    return ($ref, \@condS, \@condD);
}
sub extendLoc {
    my ($locObj, $refDb, $opt) = rearrange(['loc', 'refdb', 'opt'], @_);
    $opt = [$opt] unless ref($opt) eq "ARRAY";
    my $seqid = $locObj->seq_id;
    my $strand = $locObj->strand;
    die "not uniform seqid\n".Dumper($locObj) unless $seqid;
    die "not uniform strand\n".Dumper($locObj) unless $strand;
    my @locs = $locObj->each_Location();
    my $minS = min(map {$_->start} @locs);
    my $maxE = max(map {$_->end  } @locs);
    my $seqLen = getSeqLen($seqid, $refDb);
    for (@$opt) {
        if(/^up(\d+)$/) {
            $minS -= $1 if $strand == 1;
            $maxE += $1 if $strand == -1;
        } elsif(/^down(\d+)$/) {
            $minS -= $1 if $strand == -1;
            $maxE += $1 if $strand == 1;
        } else {
            die "unknown opt: $_\n";
        }
    }
    $minS = 1 if $minS < 1;
    $maxE = $seqLen if $maxE > $seqLen;
    my $loc = Bio::Location::Split->new();
    for ($locObj->each_Location()) {
        my ($s, $e) = ($_->start, $_->end);
        $s = $minS if $s == $locObj->min_start;
        $e = $maxE if $e == $locObj->max_end;
        $loc->add_sub_Location(Bio::Location::Simple->new(-start=>$s, -end=>$e));
    }
    $loc->seq_id($seqid);
    $loc->strand($strand);
    return $loc;
}
sub getRegionStat {
    my $self = shift;
    my ($loc, $opt) = rearrange(['loc', 'opt'], @_);
    my ($chr, $s, $e) = map {$loc->$_} qw/seq_id start end/;
    my $len = $e - $s + 1;
    my @rst;
    if($opt == 1) { #get GC counts
        my ($cntG) = (0);
        my $seqStr = $self->{biodb}->fetch_sequence(-seq_id=>$chr, -start=>$s, -end=>$e);
        while($seqStr =~ /[GC]/ig) {
            $cntG ++;
        }
        push @rst, $cntG;
    } elsif($opt == 2) { #get CDS/Intron/UTR/Intergenic counts
        my $pos = [[$s, $e]];
        my @fes = $self->getFeatures(-loc=>$loc, -types=>'mRNA');
        my (@poss, @types);
        my @typeAry = qw/CDS UTR5 UTR3 Intron Integenic/;
        print "$chr:$s..$e\n";
        my ($posMs, $posCs, $posU5s, $posU3s, $posIs, $posTs) = ([], [], [], [], [], []);
        for my $fe (@fes) {
            my ($posM, $posC, $posE, $posU5, $posU3, $posI) = getGeneStructure($fe);
            push @$posMs, @$posM;
            if(@$posC) { push @poss, @$posC; push @types, (1) x @$posC; }
            if(@$posU5) { push @poss, @$posU5; push @types, (2) x @$posU5; }
            if(@$posU3) { push @poss, @$posU3; push @types, (3) x @$posU3; }
            if(@$posI) { push @poss, @$posI; push @types, (4) x @$posI; }
        }
        if(@poss) {
            my $ref = tiling(\@poss, \@types, 1);
            for (@$ref) {
                my ($s, $e, $idx) = @$_;
                if($types[$idx] == 1) {
                    push @$posCs, [$s, $e];
                } elsif($types[$idx] == 2) {
                    push @$posU5s, [$s, $e];
                } elsif($types[$idx] == 3) {
                    push @$posU3s, [$s, $e];
                } else {
                    die "fuck me error at Localdb: 619\n" unless $types[$idx] == 4;
                    push @$posIs, [$s, $e];
                }
            }
        }
        my $cntC = @$posCs ? [posOvlp($pos, $posCs)]->[1] : 0;
        my $cntI = @$posIs ? [posOvlp($pos, $posIs)]->[1] : 0;
        my $cntU5 = @$posU5s ? [posOvlp($pos, $posU5s)]->[1] : 0;
        my $cntU3 = @$posU3s ? [posOvlp($pos, $posU3s)]->[1] : 0;
        my $cntT = @$posMs ? [posCmp($pos, $posMs)]->[4] : $len;
        die join(" ", $cntC, $cntI, $cntU5, $cntU3, $cntT, "!=", $len)." at $chr:$s-$e\n" if $cntC + $cntI + $cntU5 + $cntU3 + $cntT != $len;
        my $note = join(" ", map {$_->id."[".$_->source."][".[$_->get_tag_values("Note")]->[0]."]"} @fes);
        push @rst, ($cntC, $cntI, $cntU5, $cntU3, $cntT, $note);
    }
    push @rst, $len;
    return @rst;
}



1;
__END__
