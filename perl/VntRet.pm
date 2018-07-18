package VntRet;
use strict; use Init; use Common; use Localdb; use Readfile; use Seq;
use Path::Class; use DBI; use Data::Dumper; 
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/splitVntPerAcc/;
@EXPORT_OK = qw/recSeqByIds recSeqById recSeqByLoc 
    getVntByIds getCovByIds sumVnt
    getVntByLocs getCovByLocs/;
sub new {
    my $self = shift;
    my ($accs, $feDb, $refDb) = rearrange([qw/accs fedb refdb/], @_);
    my $ld = Localdb->new(-db=>$feDb);
    return bless {accs=>$accs, refDb=>$refDb, ld=>$ld}, 
        ref($self) || $self;
}

sub recSeqByIds {
    my $self = shift;
    my ($fi, $dir) = rearrange([qw/in out/], @_);
    my $ids = getIds($fi);
    my @types = qw/cds utr5 utr3 intron up1k dw1k/;

    my $fo = file($dir, "10_info.tbl");
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", "id", @types)."\n";
    print join("\t", "id", @types)."\n";
  
    my $dh = {};
    for my $i (0..$#types) {
        my $type = $types[$i];
        my $dirO = dir($dir, (11+$i)."_$type");
        system("mkdir -p $dirO") unless -d $dirO;
        $dh->{$type} = $dirO;
    }
  
    my $cnt = 0;
    for my $id (@$ids) {
        my $tag;
        my $seqhash = $self->recSeqById($id);
        next unless $seqhash;
        my @tags;
        for my $type (@types) {
            my $seqs = $seqhash->{$type};
            if($seqs) {
                my $fo = file($dh->{$type}, "$id.fa");
                writeSeqInOrder(-out=>$fo, -ids=>$self->{accs}, -seqs=>$seqs);
                push @tags, 1;
            } else {
                push @tags, 0;
            }
        }
        print $fh join("\t", $id, @tags)."\n";
        print join("\t", $id, @tags)."\n";
        last if ++$cnt < 1;
    }
}
sub printLoc {
    my ($seqid, $strand, $pos) = @_; 
    my ($posM, $posC, $posE, $posU5, $posU3, $posI, $posUp, $posDw) = 
        map {$pos->{$_}} qw/mrna cds exon utr5 utr3 intron up1k dw1k/;
    printf "  SeqId: %s  Strand: %d\n", $seqid, $strand;
    printf "  mRNA:\t%s\n", locAry2Str($posM);
    printf "    CDS:\t%s\n", locAry2Str($posC);
    printf "    exon:\t%s\n", locAry2Str($posE);
    printf "    utr5:\t%s\n", locAry2Str($posU5);
    printf "    utr3:\t%s\n", locAry2Str($posU3);
    printf "    intron:\t%s\n", locAry2Str($posI);
    printf "  upstream:\t%s\n", locAry2Str($posUp);
    printf "  downstream:\t%s\n", locAry2Str($posDw);
}
sub recSeqById {
    my $self = shift;
    my ($id, $opt) = @_;
    my $accs = $self->{accs};
    my @types = qw/cds utr5 utr3 intron up1k dw1k/;
    
    my $ld = $self->{ld};
    my ($seqid, $strand, $pos) = $ld->getLoc(-id=>$id, -opt=>2, -refdb=>$self->{refDb});
    my ($posM, $posC, $posU5, $posU3, $posI, $posUp, $posDw) = map {$pos->{$_}} ("mrna", @types);

    unless( $seqid =~ /^chr[1-8]$/ ) {
        print "\t$id on $seqid is skipped\n";
        return undef;
    }
    my @poss = @{$posM->[0]};
    push @poss, @{$posUp->[0]} if @$posUp;
    push @poss, @{$posDw->[0]} if @$posDw;
    my $beg = min(@poss);
    my $end = max(@poss);
    my $loc = Bio::Location::Simple->new(-start=>$beg, -end=>$end, -strand=>1, -seq_id=>$seqid);
#  printLoc($seqid, $strand, $pos);
    
    my $seqs = $self->recSeqByLoc($loc);
    my $seqhash = $self->sliceSeq($seqs, $loc, $pos);
    
    if($strand == -1) {
        for my $type (@types) {
            next unless exists $seqhash->{$type};
            my $seqs = $seqhash->{$type}; 
            for my $acc (@$accs) {
                $seqs->{$acc} = Bio::Seq->new(-seq=>$seqs->{$acc})->revcom->seq;
            }
        }
    }
    return $seqhash;
}
sub recSeqByLoc {
    my $self = shift;
    my ($loc) = @_;
    my $accs = $self->{accs};
    my $refSeq = seqRet($loc, $self->{refDb});
    my $seq_raw = { map {$_ => $refSeq} @$accs };

    die "cannot recover seq on ".$loc->seq_id."\n" if $loc->seq_id =~ /chr[^1-8]/;
    my $seq_with_cov = $self->applyCov($seq_raw, $loc);
    my $seq_with_vnt = $self->applyVnt($seq_with_cov, $loc);
    return $seq_with_vnt;
}
sub applyCov {
    my $self = shift;
    my ($seqh, $loc) = @_;
    my $accs = $self->{accs};
    my $it = getCov(-loc=>$loc, -covopt=>2, -refdb=>$self->{refDb});
    while( my $ref = $it->() ) {
        my ($chr, $pos, $posR) = map {$ref->{$_}} qw/chr pos posR/;
        my @covs = map {$ref->{cov}->{$_}} @$accs;
        for my $acc (@$accs) {
            if($ref->{cov}->{$acc} < 2) {
                substr($seqh->{$acc}, $posR-1, 1, "N");
            }
        }
    }
    return $seqh;
}
sub applyVnt {
    my $self = shift;
    my ($seqh, $loc) = @_;
    my $accs = $self->{accs};
    
    my $it = getVnt(-loc=>$loc, -refdb=>$self->{refDb});
    while( my $ref = $it->() ) {
        my ($chr, $pos, $posR, $vnt, $refa) = map {$ref->{$_}} qw/chr pos posR vnt refa/;
        for my $acc (@$accs) {
            substr($seqh->{$acc}, $posR-1, 1, $vnt->{$acc});
        }
    }
    return $seqh;
}
sub sliceSeq {
    my $self = shift;
    my ($seqs, $loc, $pos) = @_;
    my $accs = $self->{accs};
    my @types = qw/cds utr5 utr3 intron up1k dw1k/;
    my ($beg, $end) = ($loc->start, $loc->end);

    my $h;
    for my $type (@types) {
        my $posG = $pos->{$type};
        next unless @$posG;
        $posG = [ sort {$a->[0] <=> $b->[0]} @$posG ];
        my $posR = [ map {[$_->[0]-$beg+1, $_->[1]-$beg+1]} @$posG ];
        my $seqhash = { map {$_ => getSubSeq($seqs->{$_}, $posR)} @$accs };
        $h->{$type} = $seqhash;
    }
    return $h;
}

sub getVntByIds {
    my $self = shift;
    my ($fi, $fo) = rearrange([qw/in out/], @_);
    my $ids = getIds($fi);
    my $accs = $self->{accs};
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id chr pos type ref allele_cnt alleles GR MAF/, @$accs)."\n";

    my $cnt = 0;
    for my $id (@$ids) {
        $self->getVntById($id, $fh);
        last if ++$cnt < 1;
    }
    close $fh;
}
sub getVntById {
    my $self = shift;
    my ($id, $fh) = @_;
    my $accs = $self->{accs};
    my @types = qw/cds utr5 utr3 intron up1k dw1k/;
    my $ld = $self->{ld};
    my ($seqid, $strand, $pos) = $ld->getLoc(-id=>$id, -opt=>2, -refdb=>$self->{refDb});
    my ($posM, $posC, $posU5, $posU3, $posI, $posUp, $posDw) = map {$pos->{$_}} ("mrna", @types);
    
    unless( $seqid =~ /^chr[1-8]$/ ) {
        print "\t$id on $seqid is skipped\n";
        return undef;
    }
    my @poss = @{$posM->[0]};
    push @poss, @{$posUp->[0]} if @$posUp;
    push @poss, @{$posDw->[0]} if @$posDw;
    my $beg = min(@poss);
    my $end = max(@poss);
    my $loc = Bio::Location::Simple->new(-start=>$beg, -end=>$end, -strand=>1, -seq_id=>$seqid);

    my @locs;
    for my $type (@types) {
        my $loc = $pos->{$type};
        push @locs, map {[@$_, $type]} @$loc;
    }

    my $it = getVnt(-loc=>$loc, -refdb=>$self->{refDb});
    while( my $ref = $it->() ) {
        my ($chr, $pos, $posR, $vnt, $refa) = map {$ref->{$_}} qw/chr pos posR vnt refa/;

        $vnt = { map {$_=>$vnt->{$_}} @$accs };
        my @alleles = map {$vnt->{$_}} @$accs;
        my ($allele_str, $allele_num, $gr, $maf) = vntStat($vnt, $refa);
        next if $allele_num <= 1;

        my $typeIdx = first_index {$_->[0] <= $pos && $_->[1] >= $pos} @locs;
        my $type = $locs[$typeIdx]->[2];

        print $fh join("\t", $id, $chr, $pos, $type, $refa, $allele_str, $allele_num, $gr, $maf, @alleles)."\n"; 
    }
}
sub vntStat {
    my ($vnt, $refa) = @_;
    my @inds = keys %$vnt;
    my ($flag, $aCnt, $maf, $gr);
    for my $ind (@inds) {
        my $nt = $vnt->{$ind};
        if($nt !~ /^[-ATCGN]$/i) {
            die "unknown allele: $nt\n";
        } elsif($nt ne "N") {
            $aCnt->{$nt} = 0 unless exists $aCnt->{$nt};
            $aCnt->{$nt} ++;
        }
    }
    my @vars = keys %$aCnt;
    my $alleleNum = @vars;
    if($alleleNum == 0) {
        ($maf, $gr) = (0, 0);
    } else {
        $maf = $alleleNum == 1 ? 0 : min(values %$aCnt) / sum(values %$aCnt);
        $gr  = sum(values %$aCnt) / @inds;
    }
    $maf = sprintf "%.03f", $maf;
    $gr  = sprintf "%.03f", $gr;
    my $allele_str = join(" ", map {$_.":".$aCnt->{$_}} @vars);
    return ($allele_str, $alleleNum, $gr, $maf);
}
sub splitVntPerAcc {
    my ($fi, $fo) = rearrange(['in', 'out'], @_);
    my $t = readTable(-in=>$fi, -header=>1);
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id chr pos acc ref vnt/)."\n";
    my @inds = grep /^HM\d+/, $t->header;
    for my $i (0..$t->nofRow-1) {
        my ($id, $chr, $pos, $ref) = map {$t->elm($i, $_)} qw/id chr pos ref/;
        for my $ind (@inds) {
            my $var = $t->elm($i, $ind);
            if($var ne $ref && $var ne "N") {
                print $fh join("\t", $id, $chr, $pos, $ind, $ref, $var)."\n";
            }
        }
    }
}

sub getVntByLocs {
    my $self = shift;
    my ($fi, $fo) = rearrange([qw/in out/], @_);
    my $accs = $self->{accs};
    my $t = readTable(-in=>$fi, -header=>1);
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id chr pos ref allele_cnt alleles GR MAF/, @$accs)."\n";

    for my $i (0..$t->nofRow-1) {
        last if $i < 0;
        my ($id, $seqid, $beg, $end) = $t->row($i);
        my $loc = Bio::Location::Simple->new(-seq_id=>$seqid, -start=>$beg, -end=>$end, -strand=>1);
        unless( $seqid =~ /^chr[1-8]$/ ) {
            print "\t$id on $seqid is skipped\n";
            next;
        }
        $self->getVntByLoc($loc, $fh, $id);
    }
    close $fh;
}
sub getVntByLoc {
    my $self = shift;
    my ($loc, $fh, $id) = @_;
    my $accs = $self->{accs};

    my $it = getVnt(-loc=>$loc, -refdb=>$self->{refDb});
    while( my $ref = $it->() ) {
        my ($chr, $pos, $posR, $vnt, $refa) = map {$ref->{$_}} qw/chr pos posR vnt refa/;
        
        $vnt = { map {$_=>$vnt->{$_}} @$accs };
        my @alleles = map {$vnt->{$_}} @$accs;
        my ($allele_str, $allele_num, $gr, $maf) = vntStat($vnt, $refa);
        next if $allele_num <= 1;

        print $fh join("\t", $id, $chr, $pos, $refa, $allele_str, $allele_num, $gr, $maf, @alleles)."\n"; 
    }
}
sub sumVnt {
    my $self = shift;
    my ($fi, $fo, $opt) = rearrange([qw/in out opt/], @_);
    $opt ||= 1;
    my $accs = $self->{accs};
    my $t = readTable(-in=>$fi, -header=>1);
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id type acc vnt_cnt/)."\n";
    my @ids;
    if($opt == 1) {
        @ids = $t->col("id");
    } else {
        die "unknow opt $opt\n" unless $opt == 2;
        @ids = map {join("+", $t->elm($_, "id"), $t->elm($_, "type"))} (0..$t->nofRow-1);
    } 
    my $ref = group(\@ids);
    for my $id (sort keys %$ref) {
        my ($idx, $nofId) = @{$ref->{$id}};
        my @idxs = ($idx..$idx+$nofId-1);
        my $snph = { map {$_=>0} @$accs };
        for my $i ($idx..$idx+$nofId-1) {
            my $refa = $t->elm($i, "ref");
            for my $acc (@$accs) {
                my $alta = $t->elm($i, $acc);
                $snph->{$acc} ++ if $alta ne "N" && $alta ne $refa;
            }
        }
        my $type = "";
        ($id, $type) = split(/\+/, $id) if $opt == 2;
        for my $acc (@$accs) {
            print $fh join("\t", $id, $type, $acc, $snph->{$acc})."\n";
        }
    }
    close $fh;
}

sub getCovByIds {
    my $self = shift;
    my ($fi, $fo) = rearrange([qw/in out/], @_);
    my $ids = getIds($fi);
    my $accs = $self->{accs};
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id type chr loc length acc nN cov1_mean cov1_median cov1_sd cov2_mean cov2_median cov2_sd/)."\n";

    my $cnt = 0;
    for my $id (@$ids) {
        print "  $id YES\n" if $self->getCovById($id, $fh);
        last if ++$cnt == 0;
    }
    close $fh;
}
sub getCovById {
    my $self = shift;
    my ($id, $fh) = @_;
    my $accs = $self->{accs};
    my @types = qw/cds utr5 utr3 intron up1k dw1k/;
  
    my $ld = $self->{ld};
    my ($seqid, $strand, $pos) = $ld->getLoc(-id=>$id, -opt=>2, -refdb=>$self->{refDb});
    my ($posM, $posC, $posU5, $posU3, $posI, $posUp, $posDw) = map {$pos->{$_}} ("mrna", @types);
    
    unless( $seqid =~ /^chr[1-8]$/ ) {
        print "\t$id on $seqid is skipped\n";
        return undef;
    }
    my @poss = @{$posM->[0]};
    push @poss, @{$posUp->[0]} if @$posUp;
    push @poss, @{$posDw->[0]} if @$posDw;
    my $beg = min(@poss);
    my $end = max(@poss);
    my $loc = Bio::Location::Simple->new(-start=>$beg, -end=>$end, -strand=>1, -seq_id=>$seqid);
    
    my $covh1 = getRegionCov($loc, 1, $accs, $self->{refDb});
    my $covh2 = getRegionCov($loc, 2, $accs, $self->{refDb});
    my $co_cov2 = 2;
    for my $acc (@$accs) {
        my ($cov1, $cov2) = map {$_->{$acc}} ($covh1, $covh2);
        for my $type (@types) {
            my $locAry = $pos->{$type};
            next unless @$locAry;
            my $len = locAryLen($locAry);
            my ($locStr) = locObj2Str(locAry2Obj($locAry, $strand, $seqid));

            my $locAryR = [ map { [$_->[0]-$beg+1, $_->[1]-$beg+1] } @$locAry ];
            my $cov_s1 = [ map {@$cov1[$_->[0]-1..$_->[1]-1]} @$locAryR ];
            my $cov_s2 = [ map {@$cov2[$_->[0]-1..$_->[1]-1]} @$locAryR ];
            die "cov_s1 [".@$cov_s1."] not $len bp\n" unless @$cov_s1 == $len;
            die "cov_s2 [".@$cov_s2."] not $len bp\n" unless @$cov_s2 == $len;
            my @stats1 = getDescStat($cov_s1);
            my @stats2 = getDescStat($cov_s2);
            my $nN = grep {$_ < $co_cov2} @$cov_s2;
            print $fh join("\t", $id, $type, $seqid, $locStr, $len, $acc, $nN, @stats1, @stats2)."\n";
        }
    }
    return 1;
}

sub getCovByLocs {
    my $self = shift;
    my ($fi, $fo) = rearrange([qw/in out/], @_);
    my $accs = $self->{accs};
    my $t = readTable(-in=>$fi, -header=>1);
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id chr loc length acc nN cov1_mean cov1_median cov1_sd cov2_mean cov2_median cov2_sd/)."\n";

    for my $i (0..$t->nofRow-1) {
        last if $i < 0;
        my ($id, $seqid, $beg, $end) = $t->row($i);
        my $loc = Bio::Location::Simple->new(-seq_id=>$seqid, -start=>$beg, -end=>$end, -strand=>1);
        my ($locStr) = locObj2Str($loc);
        unless( $seqid =~ /^chr[1-8]$/ ) {
            print "\t$id on $seqid is skipped\n";
            next;
        }
        my $len = locObjLen($loc);

        my $covh1 = getRegionCov($loc, 1, $accs, $self->{refDb});
        my $covh2 = getRegionCov($loc, 2, $accs, $self->{refDb});
        my $co_cov2 = 2;
        for my $acc (@$accs) {
            my ($cov1, $cov2) = map {$_->{$acc}} ($covh1, $covh2);
            my @stats1 = getDescStat($cov1);
            my @stats2 = getDescStat($cov2);
            my $nN = grep {$_ < $co_cov2} @$cov2;
            print $fh join("\t", $id, $seqid, $locStr, $len, $acc, $nN, @stats1, @stats2)."\n";
        }
    }
    close $fh;
}


1;
__END__
