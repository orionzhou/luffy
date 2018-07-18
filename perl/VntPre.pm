package Vnt;
use strict; use Init; use Common; use Localdb; use Readfile;
use Path::Class; use DBI; use IO::File; use Data::Dumper; use Parser;
use Clone qw/clone/; use Writefile; use DB_File;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Math::Round qw/round/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/reformatIndel1 reformatIndel2 parse_vnt_line parse_vnt_file
    recover_vnt_line getAccFromHead storeVnt vntDbHead vntDbLoop vntDbPos getVntStat
    vntFilter1 vntFilter2 callVnt storeCalledVnt getCalledVnt/;
@EXPORT_OK = qw//;
sub parse_vnt_line {
    my ($line, $param) = @_;
    my $paramD = {colI=>0, colA=>6, lenA=>6, cntA=>31};
    my ($colI, $colA, $lenA, $cntA) = 
        map {exists $param->{$_} ? $param->{$_} : $paramD->{$_}} qw/colI colA lenA cntA/;
    my @ps = split("\t", $line);
    die "< not enough cols at\n $line\n" if @ps < $colA+$lenA*$cntA;
    my $ref = {
        seqid => $ps[$colI+0],
        id    => $ps[$colI+1],
        class => $ps[$colI+2],
        pos2  => $ps[$colI+3],
        ref   => $ps[$colI+4],
        var   => $ps[$colI+5],
        line  => $line
    };
    my $pos;
    my $class = $ref->{class};
    if($class eq "S") {
        $pos = $1 if $ref->{pos2} =~ /^(\d+)$/;
    } elsif($class eq "I") {
        $pos = $1 if $ref->{pos2} =~ /^(\d+)\^(\d+)$/;
        $pos = $1 if $ref->{pos2} =~ /^([\d\.]+)$/;
    } elsif($class eq "D") {
        $pos = $1 if $ref->{pos2} =~ /^\[(\d+)\.\.(\d+)\]$/;
        $pos = $1 if $ref->{pos2} =~ /^(\d+)$/;
    } else {
        die "Unknown class [$class] at ".$ref->{pos2}."\n";
    }
    die "illegal position format: ".$ref->{pos2}."\n" unless $pos;
    $ref->{pos} = $pos;
    my $refa;
    for my $i (0..$cntA-1) {
        my $colAcc = $colA + $lenA * $i;
        my $acc = sprintf("HM%03d", $i+1);
        $acc = "HM101" if $i == 30;
        $refa->{$acc} = {
            reads1   => $ps[$colAcc+0] ? $ps[$colAcc+0] : 0,
            cov1     => $ps[$colAcc+1] ? $ps[$colAcc+1] : 0,
            reads2   => $ps[$colAcc+2] ? $ps[$colAcc+2] : 0,
            cov2     => $ps[$colAcc+3] ? $ps[$colAcc+3] : $ps[$colAcc+2],
            avgQual  => $ps[$colAcc+4] ? $ps[$colAcc+4] : 0,
            maxQual  => $ps[$colAcc+5] ? $ps[$colAcc+5] : 0};
    }
    $ref->{acc} = $refa;
    return $ref;
}
sub parse_vnt_file {
    my ($fi, $accs, $param) = rearrange(["in", "acc", "param"], @_);
    die "$fi is not there\n" unless -s $fi;
    my $fh = new IO::File $fi, "r";
    return sub {
        if( eof($fh) ) {
            return undef;
        } else {
            my $line = "";
            while(!$line) {
                $line = readline($fh);
                $line =~ s/\r?\n?$//;
                $line = "" if $line =~ /^Reference/i;
            }
            return parse_vnt_line($line, $param);
        }
    }
}
sub recover_vnt_line {
    my ($ref) = @_;
    my $refa = $ref->{acc};
    my @ps = map {$ref->{$_}} qw/seqid id class pos2 ref var/;
    for my $acc (sort keys %$refa) {
        my $refacc = $refa->{$acc};
        push @ps, map {$refacc->{$_}} qw/reads1 cov1 reads2 cov2 avgQual maxQual/;
    }
    return join("\t", @ps);
}
sub getAccFromHead {
    my ($fi) = @_;
    my $fh = new IO::File $fi, "r";
    my $tmp = readline($fh);
    chomp $tmp;
    my @ps = split("\t", $tmp);
    close $fh;
    my $ref;
    for my $i (0..$#ps) {
        if($ps[$i] =~ /^(HM[0-9]+)/i) {
            $ref->{$1} ||= $i;
        }
    }
    return $ref;
}
sub reformatIndel1 {
#deal with multiple basepair insertion polymorphism
#change "100^103 I --- ATT" to "[101..102] D GT --" and "100^101 I - ATT"
    my ($seqid, $refDb) = @_;
    my $fi = file($DIR_Variant, "01_raw", "$seqid.txt");
    my $fo = file($DIR_Variant, "02_reformat", "$seqid.txt");
    die "$fi is not there\n" unless -s $fi;
    my $ld = Localdb->new(-db=>$refDb);
    my $fh = new IO::File $fo, "w";
    my $posPrev = 0;
    my $ps = parse_vnt_file(-in=>$fi, -param=>{colA=>16});
    while(my $ref = $ps->()) {
        my ($seqid, $pos, $pos2, $class) = map {$ref->{$_}} qw/seqid pos pos2 class/;
        if($seqid =~ /chr(\d+)/i) {
            $ref->{seqid} = "chr$1";
        } elsif($seqid =~ /mitochon/i) {
            $ref->{seqid} = "mitochondrion";
        } elsif($seqid =~ /AC093544/i) {
            $ref->{seqid} = "chloroplast";
        } elsif($seqid =~ /^([A-Za-z]{1,3}\_?[\d\.]{3,10})$/) {
            my $pos2 = $ref->{pos2};
            if($pos2 =~ /^\[?(\d+)(\^|\.\.)(\d+)\]?$/) {
                my $loc1 = "$seqid:$1..$3";
                my $loc2 = $ld->bac2Chr($loc1);
                my ($seqid2, $start2, $end2) = splitLocStr($loc2);
                $ref->{seqid} = $seqid2;
                $ref->{pos2} = join($2, $start2, $end2);
                $ref->{pos2} = "[".$ref->{pos2}."]" if $2 eq "..";
            } else {
                die "unknown position flag: $pos2\n" unless $pos2 =~ /^\d+$/;
                my $loc2 = $ld->bac2Chr("$seqid:$pos2..$pos2");
                my ($seqid2, $start2, $end2) = splitLocStr($loc2);
                $ref->{seqid} = $seqid2;
                $ref->{pos2} = $start2;
            }
        }
        $ref->{id} = $1 if $ref->{id} =~ /(\d+)/;
        my $flag = 0;
        if($class eq "I") {
            $pos2 =~ /^(\d+)\^(\d+)$/;
            my ($start, $end) = ($1, $2);
            if($start+1 < $end) {
                $flag = 1;
                my $var = $ref->{var};
                #deletion
                my $locStr = $ref->{seqid}.":".($start+1)."..".($end-1);
                my $seqStr = seqRet($locStr, $refDb);
                $ref->{class} = "D";
                $ref->{pos2}  = "[".($start+1)."..".($end-1)."]";
                $ref->{ref}   = $seqStr;
                $ref->{var}   = "-" x ($end-$start-1);
                print $fh recover_vnt_line($ref)."\n";
                #insertion
                $ref->{class} = "I";
                $ref->{pos2}  = $start."^".($start+1);
                $ref->{ref}   = "-" x length($var);
                $ref->{var}   = $var;
                print $fh recover_vnt_line($ref)."\n";
            } else {
                die "illegal indel: start[$start] - end[$end]\n" unless $start+1 == $end;
            }
        }
        print $fh recover_vnt_line($ref)."\n" if $flag == 0;
    }
}
sub reformatIndel2 {
    my ($seqid) = @_;
    my $fi = file($DIR_Variant, "02_reformat", "$seqid.txt");
    my $fo = file($DIR_Variant, "03_reformat", "$seqid.txt");
    die "$fi is not there\n" unless -s $fi;
    my $fh = new IO::File $fo, "w";
    my $ps = parse_vnt_file(-in=>$fi);
    my $cntVnt = 0;
    my $cache = {};
    while(my $ref = $ps->()) {
        my ($seqid, $pos, $pos2, $class) = map {$ref->{$_}} qw/seqid pos pos2 class/;
        my $ref2 = clone($ref);
        die "unsupported class: $class at $pos2\n" unless $class =~ /[SDI]/;
        if($class eq "S") {
            $cache->{$pos} ||= [];
            my $ref3s = $cache->{$pos};
                my $idx = first_index {$ref3s->[$_]->{var} eq $ref2->{var}} (0..@$ref3s-1);
                die "$seqid:$pos has 2 S allels\n" if $idx != -1;
                push @$ref3s, $ref2;
        } elsif($class eq "D") {
            $pos2 =~ /^\[(\d+)\.\.(\d+)\]$/;
            my ($start, $end, $delLen) = ($1, $2, $2-$1+1);
            my $refLen = length($ref->{ref});
            die "refLen[$refLen] != delLen[$delLen] at $pos2\n" unless $delLen == $refLen;
            my $step = $delLen<10 ? 0.1 : 0.01;
            for my $i (0..$delLen-1) {
                my $ref2 = clone($ref);
                my $posD = $pos + $i;
                $ref2->{pos2} = $posD;
                $ref2->{pos}  = $posD;
                $ref2->{ref}  = substr($ref->{ref}, $i, 1);
                $ref2->{var}  = '-';
                $ref2->{id}   = $ref->{id} + $i*$step;
                $cache->{$posD} ||= [];
                my $ref3s = $cache->{$posD};
                my $idx = first_index {$ref3s->[$_]->{var} eq "-"} (0..@$ref3s-1);
                if($idx == -1) {
                    push @$ref3s, $ref2;
                } else {
                    my $ref3 = $ref3s->[$idx];
                    $ref3 = infoUpate($ref3, $ref2);
                }
            }
        } else {
            my $insLen = length($ref->{var});
            my $refLen = length($ref->{ref});
            die "refLen[$refLen] != insLen[$insLen] at $pos2\n" unless $insLen == $refLen;
            my $step = $insLen<10 ? 0.1 : 0.01;
            for my $i (1..$insLen) {
                my $ref2 = clone($ref);
                my $posI = $pos + $i*$step;
                $ref2->{pos2} = $posI;
                $ref2->{pos} = $posI;
                $ref2->{ref} = "-";
                $ref2->{var} = substr($ref->{var}, $i-1, 1);
                $ref2->{id} = $ref->{id} + $i*$step;
                $cache->{$posI} ||= [];
                my $ref3s = $cache->{$posI};
                my $idx = first_index {$ref3s->[$_]->{var} eq $ref2->{var}} (0..@$ref3s-1);
                if($idx == -1) {
                    push @$ref3s, $ref2;
                } else {
                    my $ref3 = $ref3s->[$idx];
                    $ref3 = infoUpate($ref3, $ref2);
                }
            }
        }
        $cache = dumpCache($cache, $fh);
        last if ++$cntVnt < 0;
    }
    $cache = dumpCache($cache, $fh, "all");
    if($seqid eq 'chr0') {
        print "\tsorting file using chr position for $seqid...\n";
        my $fTmp = $fo.".tmp";
        sortFileByCol(-in=>$fo, -out=>$fTmp, -cols=>'4n');
        system("mv $fTmp $fo");
    }
    print "$seqid: $cntVnt variants processed\n";
}
sub infoUpate {
    my ($ref3, $ref4) = @_;
    my ($ref1, $ref2) = (clone($ref3), clone($ref4));
    my ($refO, $refN) = map {$_->{ref}} ($ref1, $ref2);
    my ($varO, $varN) = map {$_->{var}} ($ref1, $ref2);
    die "oldRef[$refO] != newRef[$refN]\n".$ref1->{line}."\n".$ref2->{line}."\n" unless $refO eq $refN;
    die "oldVar[$varO] != newVar[$varN]\n".$ref1->{line}."\n".$ref2->{line}."\n" unless $varO eq $varN;
    $ref1->{id} .= "+".$ref2->{id} if $ref1->{Id} ne $ref2->{id};
    my ($refAcc1, $refAcc2) = map {$_->{acc}} ($ref1, $ref2);
    for my $acc (keys %$refAcc1) {
        my ($ref5, $ref6) = map {$_->{$acc}} ($refAcc1, $refAcc2);
        my ($numOld, $numNew) = ($ref5->{reads1}, $ref6->{reads1});
        $ref5->{cov1} = max($ref5->{cov1}, $ref6->{cov1});
        $ref5->{cov2} = max($ref5->{cov2}, $ref6->{cov2});
        $ref5->{reads1} += $ref6->{reads1};
        $ref5->{reads2} += $ref6->{reads2};
        if($ref5->{reads1} > 0) {
            my $avgOld = $ref5->{avgQual} ? $ref5->{avgQual} : 0;
            my $avgNew = $ref6->{avgQual} ? $ref6->{avgQual} : 0;
            my $maxOld = $ref5->{maxQual} ? $ref6->{maxQual} : 0;
            my $maxNew = $ref5->{maxQual} ? $ref6->{maxQual} : 0;
            $ref5->{maxQual} = max($maxOld, $maxNew);
            $ref5->{avgQual} = sprintf("%.02f", ($avgOld*$numOld + $avgNew*$numNew) / ($numOld+$numNew) );
        }
    }
    return $ref1;
}
sub dumpCache {
    my ($cache, $fh, $option) = @_;
    my $dumpHash = {};
    my $cacheSize = 200;
    my $cacheSizeNow = scalar(keys(%$cache));
    if($option eq "all") {
        $dumpHash = $cache;
        die "cacheSize[$cacheSizeNow] not $cacheSize" unless $cacheSizeNow == $cacheSize;
        $cache = {};
    } else {
        my $cnt = $cacheSizeNow - $cacheSize;
        if($cnt>0) {
            #die "not 1 more : $cnt\n" unless $cnt == 1;
            my @sortedKeys = sort({$a<=>$b} keys(%$cache));
            for(my $i=0; $i<$cnt; $i++) {
                $dumpHash->{$sortedKeys[$i]} = $cache->{$sortedKeys[$i]};
                delete $cache->{$sortedKeys[$i]};
            }
        }
    }
    for my $pos (sort({$a<=>$b} keys(%$dumpHash))) {
        my $refs = $dumpHash->{$pos};
        for my $ref (@$refs) {
            print $fh recover_vnt_line($ref)."\n";
        }
    }
    return $cache;
}
sub numCompare {
    my ($k1, $k2) = @_;
    return $k1 <=> $k2;
}
sub storeVnt {
    my ($seqid) = @_;
    my $fi = file($DIR_Variant, "03_reformat", "$seqid.txt");
    my $fo = file($DIR_Variant, "04_vnt_db", "$seqid.db");
    my %h;
    unlink $fo;
    $DB_BTREE->{'compare'} = \&numCompare;
    $DB_BTREE->{'flags'} = R_DUP;
    tie %h, "DB_File", $fo, O_RDWR|O_CREAT, 0666, $DB_BTREE
        or die "cannot bind to $fo\n:$!\n";
    my $ps = parse_vnt_file($fi);
    while(my $ref = $ps->()) {
        my $key = $ref->{pos};
        $h{$key} = $ref->{line};
    }
    untie %h;
}
sub vntDbHead {
    my ($sid) = @_;
    my $f = file($DIR_Variant, "04_vnt_db", "$sid.db");
    my %h;
    $DB_BTREE->{'compare'} = \&numCompare;
    $DB_BTREE->{'flags'} = R_DUP;
    my $db = tie %h, "DB_File", $f, O_RDONLY, 0444, $DB_BTREE
        or die "cannot bind to $f\n:$!\n";
    my ($k, $v);
    my $cnt = 0;
    my $st = $db->seq($k, $v, R_FIRST);
    while($st == 0 && ++$cnt < 20) {
        my $ref = parse_vnt_line($v);
        print "\t".join("\t", $k, $ref->{ref}, $ref->{class}, $ref->{var})."\n";
        $st = $db->seq($k, $v, R_NEXT);
    }
    untie %h;
}
sub vntDbLoop {
    my ($loc) = @_;
    my ($sid, $start, $end) = splitLocStr($loc);
    my $f = file($DIR_Variant, "04_vnt_db", "$sid.db");
    my %h;
    $DB_BTREE->{'compare'} = \&numCompare;
    $DB_BTREE->{'flags'} = R_DUP;
    my $db = tie %h, "DB_File", $f, O_RDONLY, 0444, $DB_BTREE
        or die "cannot bind to $f\n:$!\n";
    my ($k, $v) = ($start, 0);
    my ($st, $cnt) = (0, 0);
    return sub {
        if( ++$cnt == 1 ) {
            $st = $db->seq($k, $v, R_CURSOR);
        } else {
            $st = $db->seq($k, $v, R_NEXT);
        }
        if($st == 0 && $k <= $end) {
            return [$k, parse_vnt_line($v)];
        } else {
            return undef;
        }
    }
}
sub vntDbPos {
    my ($loc) = @_;
    my $l = vntDbLoop($loc);
    my $rst;
    while(my $r = $l->() ) {
        my ($pos, $ref) = @$r;
        $rst->{$pos} ||= [];
        delete $ref->{line};
        push @{$rst->{$pos}}, $ref;
    }
    return $rst;
}
sub getVntStat {
    my ($refs) = @_;
    my $rst;
    my @sids = uniq( map {$_->{seqid}} @$refs );
    die ">1 seqids: ".join(" ", @sids)."\n" unless @sids == 1;
    my @refas = uniq( map {$_->{ref}} @$refs );
    die ">1 references: ".join(" ", @refas)."\n" unless @refas == 1;
    my @poss = uniq( map {$_->{pos}} @$refs );
    die ">1 poss: ".join(" ", @poss)."\n" unless @poss == 1;
    my @accs = sort keys %{$refs->[0]->{acc}};
    for my $acc (@accs) {
        for my $ref (@$refs) {
            my $stat = {
                var    => $ref->{var},
                class  => $ref->{class}, 
                cov1   => $ref->{acc}->{$acc}->{cov1},
                cov2   => $ref->{acc}->{$acc}->{cov2},
                reads1 => $ref->{acc}->{$acc}->{reads1},
                reads2 => $ref->{acc}->{$acc}->{reads2}
            };
            $stat->{freq} = $stat->{cov1} > 0 ? sprintf( "%.03f", $stat->{reads1} / $stat->{cov1} ) : 0;
            $rst->{$acc} ||= [];
            push @{$rst->{$acc}}, $stat;
        }
    }
    return ($rst, $refas[0]);
}
sub vntFilter1 {
    my ($ref, $co) = @_;
    $co ||= {freq=>0.7, reads2=>2};
    my $rst;
    for my $acc (keys %$ref) {
        my $ref2 = $ref->{$acc};
        push @$ref2, {
            freq   => 1 - sum( map {$_->{freq}} @$ref2 ),
            reads2 => max( map {$_->{cov2}} @$ref2 ) - sum( map {$_->{reads2}} @$ref2 )
        };
        my $idx = first_index {$_->{freq} >= $co->{freq} && $_->{reads2} >= $co->{reads2}} @$ref2;
        if($idx == -1) {
            $rst->{$acc} = "N";
        } elsif($idx == @$ref2-1) {
            $rst->{$acc} = "=";
        } else {
            my $var = $ref2->[$idx]->{var};
            $rst->{$acc} = $var =~ /[ATCG\-]/i ? $var : "N";
        }
    }
    return $rst;
}
sub vntFilter2 {
    my ($ref) = @_;
    my $accs = [ keys %$ref ];
    my ($rst, $flag, $aCnt, $maf, $gt);
    for my $acc (keys %$ref) {
        my $nt = $ref->{$acc};
        if($nt ne "N") {
            $aCnt->{$nt} = 0 unless exists $aCnt->{$nt};
            $aCnt->{$nt} ++;
        }
    }
    my @vars = grep /[-atgc]/i, keys %$aCnt;
    my $alleleNum = scalar(keys %$aCnt);
    if($alleleNum == 0) {
        $flag = "mono1";
        ($maf, $gt) = (0, 0);
    } else {
        $maf = min(values %$aCnt) / sum(values %$aCnt);
        $gt  = sum(values %$aCnt) / @$accs;
        if($alleleNum == 1) {
            $flag = @vars == 0 ? "mono1" : "mono2";
        } elsif($alleleNum == 2) {
            $flag = "bi";
        } elsif($alleleNum > 2) {
            $flag = "multi";
        }
    }
    $rst = {maf=>sprintf("%.03f", $maf), gt=>sprintf("%.03f", $gt), vars=>\@vars, flag=>$flag};
    return $rst;
}
sub callVnt {
    my ($seqid, $co, $refDb) = rearrange(['seqid', 'co', 'refdb'], @_);
    my $fi = file($DIR_Variant, "04_vnt_db", "$seqid.db");
    my $fo = file($DIR_Variant, "05_vntcalled", "$seqid.txt");
    my $fh = new IO::File $fo, "w";
    my $win = getWindows(-chr=>$seqid, -db=>$refDb, -winsize=>10_000, -winstep=>10_000, -opt=>2);
    my $cntVnt = 0;
    for my $i (0..@$win-1) {
        my ($wS, $wE) = @{$win->[$i]};
        my $loc = "$seqid:$wS..$wE";
        my $vnts = vntDbPos($loc);
        for my $pos (sort {$a<=>$b} keys %$vnts) {
            my $refs = $vnts->{$pos};
            my ($ref2, $refa) = getVntStat($refs);
            my $ref3 = vntFilter1($ref2, $co);
            my @ary = map {$ref3->{$_}} sort keys %$ref3;
            my $ref4 = vntFilter2($ref3, $refa);
            if($pos == 0) {
                die join("\t", $pos, $refa, prettyStr(join("", @ary)))."\n";
            }
            if($refa =~ /[ATCG\-]/i && $ref4->{flag} ne "mono1") {
                print $fh join("\t", $pos, $refa, prettyStr(join("", @ary)))."\n";
                $cntVnt ++;
                die if $cntVnt < 1;
            }
        }
    }
    print join("\t", $seqid, $cntVnt)."\n";
}
sub storeCalledVnt {
    my ($seqid) = @_;
    my $fi = file($DIR_Variant, "05_vntcalled", "$seqid.txt");
    my $fo = file($DIR_Variant, "06_vntcalled_db", "$seqid.db");
    my %h;
    unlink $fo;
    $DB_BTREE->{'compare'} = \&numCompare;
    $DB_BTREE->{'flags'} = R_DUP;
    tie %h, "DB_File", $fo, O_RDWR|O_CREAT, 0666, $DB_BTREE
        or die "cannot bind to $fo\n:$!\n";
    my $fh = new IO::File $fi, "r";
    while( <$fh> ) {
        chomp;
        my @ps = split("\t", $_, 2);
        my $k = $ps[0];
        my $v = $ps[1];
        $v =~ s/ //g;
        $h{$k} = $v;
    }
    untie %h;
}
sub getCalledVnt {
    my ($loc) = @_;
    my $sid = $loc->seq_id;
    die "not uniform seqid\n".Dumper($loc) unless $sid;
    my $f = file($DIR_Variant, "06_vntcalled_db", "$sid.db");
    my %h;
    $DB_BTREE->{'compare'} = \&numCompare;
    $DB_BTREE->{'flags'} = R_DUP;
    my $db = tie %h, "DB_File", $f, O_RDONLY, 0444, $DB_BTREE
        or die "cannot bind to $f\n:$!\n";
    my @locs = $loc->each_Location();
    my $minS = min(map {$_->start} @locs);
    my $maxE = max(map {$_->end  } @locs);
    my ($k, $v) = ($minS, 0);
    my ($st, $cnt) = (0, 0);
    return sub {
        if( ++$cnt == 1 ) {
            $st = $db->seq($k, $v, R_CURSOR);
        } else {
            $st = $db->seq($k, $v, R_NEXT);
        }
        if($st == 0 && $k <= $maxE) {
            my $idx = first_index {$k >= $_->start && $k <= $_->end} @locs;
            return 2 if $idx == -1;
            my ($ref, $var) = split("\t", $v);
            my @vars = split("", $var);
            die "not 31 accs: $k - $v\n" unless @vars == 31;
            my %vh = map {($_==30 ? "HM101" : sprintf("HM%03d", $_+1)) => $vars[$_]} (0..30);
            $vh{HM101} = "N" if $vh{HM101} ne "=";
            return [$k, $ref, \%vh];
        } else {
            return undef;
        }
    }
}

1;
__END__
