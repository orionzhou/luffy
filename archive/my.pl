#!/usr/bin/perl -w
use Bio::Seq; use Bio::SeqIO; use Bio::SearchIO; 
use Bio::DB::BioDB; use Bio::DB::GFF; use Bio::DB::SeqFeature::Store;
use Bio::Search::Tiling::MapTiling;
use Path::Class; use DBI; use Carp::Assert; use FileHandle;
use strict;

sub constants_init() {
    our $DIR_HOME;
    if($^O =~ /win32/i) {
        $DIR_HOME = dir($ENV{'HOMEPATH'});
    } else {
        our $DIR_site = dir($ENV{'HOME'}, "Sites");
        if($ENV{"mode"} == 2) {
            $DIR_HOME = dir($ENV{'HOME'});
        } else {
            $DIR_HOME = dir($ENV{'WORK'});
        }
    }
    our $DIR_data = dir($DIR_HOME, "Data");
    our $DIR_input = dir($DIR_data, "in");
    our $DIR_output = dir($DIR_data, "out");
    our $DIR_stat = dir($DIR_data, "stat");
}
sub db_init() {
    my $DSN = 'dbi:mysql:database='.$ENV{"MYSQL_DB"}.';host='.$ENV{"MYSQL_HOST"};
    my $ADAPTOR     = 'DBI::mysql';
    my $USER     = $ENV{'MYSQL_USER'};
    my $PASS     = $ENV{'MYSQL_PW'};
    our $DB = Bio::DB::SeqFeature::Store->new(
        -adaptor => $ADAPTOR, 
        -dsn     => $DSN,
        -user    => $USER, 
          -pass    => $PASS)
            or die "Couldn't create connection to database";
    our $DB_MYSQL = DBI->connect(
        "DBI:mysql:$ENV{'MYSQL_DB'}:$ENV{'MYSQL_HOST'}",
        $ENV{"MYSQL_USER"}, 
        $ENV{"MYSQL_PW"}, 
        {"RaiseError"=>1});
}
sub seq_ret_old() {
    my $DIR_data = $_;
    my ($chr, $chr_start, $chr_stop, $flkg_len) = @_;
    $flkg_len=0 if !$flkg_len;
    if( $chr_start !~ /^[0-9]{1,10}$/ || $chr_stop !~ /^[0-9]{1,10}$/  || $chr_start>=$chr_stop ) {
        die("Invalid input parameter!");
    } else {
        my $fSeq = file($DIR_data, "Genome_Mt", "MtChr".$chr);
        my $chrFH = new Bio::SeqIO(-format => 'fasta', -file => $fSeq);
        my $chrSeq = $chrFH->next_seq();
        my $length = $chr_stop - $chr_start + 1;
        my $seq_got = $chrSeq->subseq($chr_start-$flkg_len,$chr_stop+$flkg_len);
        return $seq_got;
    }
}
sub getSeq {
    my ($id, $type) = @_;
    my @features = $DB->get_features_by_name($id);
    my $seqId = join("_", $id, $type);
    my $seqStr = "";
    if(@features == 0) {
        print "No features named $id\n";
        exit 1;
    } elsif(@features > 1) {
        print "more than 1 features named $id\n";
        exit 1;
    } else{
        my $fe = $features[0];
        print $fe;
        if($fe->primary_tag eq "mRNA" && $type eq "CDS") {
            my @cdsAry = $fe->get_SeqFeatures("CDS");
            my %cdsHash;
            for my $cds (@cdsAry) {
                $cdsHash{$cds->start} = $cds;
            }
            my @cdsStartAry;
            if ($fe->strand == 1) {
                @cdsStartAry = sort keys %cdsHash;
            } else {
                assert($fe->strand == -1);
                @cdsStartAry = reverse sort keys %cdsHash;
            }
            for my $start (@cdsStartAry) {
                $seqStr .= $cdsHash{$start}->seq->seq;
            }
        } else {
            $seqStr = $fe->seq->seq;
        }
    }
    return Bio::Seq->new(-id=>$seqId, -seq=>$seqStr);
}
sub getLocIndex_old {
    my ($id, $pos) = rearrange(["id", "pos"], @_);
    my $nameToAcc = {"MtChloro" => 'AC093544', "MtMito" => 'Y08501'};
    $id = exists($nameToAcc->{$id}) ? $nameToAcc->{$id} : $id;
    my ($query, $sth);
    my $filePos = -1;
    if($id =~ /^((MtChr\d)|[A-Za-z]{1,3}[\d\.]{3,9})$/) {
        $pos = $pos ? $pos : 1;
        if($id =~ /^([A-Za-z]{1,3}\d{3,6})\.(\d+)$/) {
            $id = $1;
        }
        $query = "SELECT id FROM mt_hapmap.mt32 WHERE seqid=? and start<=? and type='mRNA' ORDER BY start DESC LIMIT 1";
        $sth = $dbh_mthap->prepare($query);
        $sth->execute($id, $pos) or die "Couldn't execute query: " . $sth->errstr;
        if($sth->rows == 1) {
            my @tmp = $sth->fetchrow_array();
            $id = $tmp[0];
        } elsif($sth->rows == 0) {
            return $filePos;
        }
    }
    $query = 'SELECT note FROM mt_hapmap.mt32 WHERE id = ?';
    $sth = $dbh_mthap->prepare($query);
    $sth->execute($id) or die "Couldn't execute query: " . $sth->errstr;
    if($sth->rows != 1) {
        print $id." not found in mt_hapmap.mt32\n";
        exit 1;
    }
    my @rsAry = $sth->fetchrow_array();
    my ($note) = ($rsAry[0]);
    $filePos = $1 if $note =~ /filePos=(\d+)/i;
    return $filePos;
}
=usage
my ($covHash1, $covHash2) = getCovS(-acc=>["HM001"], -chr=>1, -start=>100000, -end=>100100, -covopt=>1, -covthresh=>2);
=cut
sub getCovS {
    my ($accAry, $chr, $start, $end, $covOpt, $cutOffCov) =
        rearrange(["acc", "chr", "start", "end", "covopt", "covthresh"], @_);
    my $chrName = $chr=~/^\d$/ ? "MtChr".$chr : $chr;
    my $totalLen = $end - $start + 1;
    $covOpt = $covOpt ? $covOpt : 1;
    $cutOffCov = $cutOffCov ? $cutOffCov : 2;
    my $covOptHash = {1=>'cov', 2=>'uniqCov'};
    my $covOption = $covOptHash->{$covOpt};
    $accAry = $accAry ? $accAry : getAccInDir(dir($DIR_Mt, "Coverage"));
    $accAry = ref($accAry) ? $accAry : [$accAry];
    my ($covHash1, $covHash2);
    my $qryExp = join("_", "^(".join("|", @$accAry).")", $chrName, $covOption);
    
    my ($wStart, $wEnd) = ( ceil( $start / 2000 ), ceil( $end / 2000 ) );
    my $offSetStart = ($wStart - 1) * 2000;
    my $offSetEnd   = ($wEnd   - 1) * 2000;
    
    my $qry = qq{SELECT A.name, B.offset, B.cov FROM mt_hapmap.cov_name AS A, mt_hapmap.cov AS B
        WHERE A.id = B.id AND A.name REGEXP ? AND B.offset >= ? AND B.offset <= ?};
    my $sth = $dbh_mthap->prepare($qry);
    $sth->execute($qryExp, $offSetStart, $offSetEnd);
    
    my $rs = {};
    while( my ($name, $offset, $cov) = $sth->fetchrow_array() ) {
        my @tmp = split("_", $name);
        $rs->{$tmp[0]}->{$offset} = $cov;
    }
    for my $acc (sort keys(%$rs)) {  
        my $covStr = "";
        for my $offset (sort {$a<=>$b} keys %{$rs->{$acc}}) {
            my $subStart = $offset+1<$start  ? $start-$offset : $offset-$offset+1;
            my $subEnd   = $offset+2000>$end ? $end-$offset   : 2000;
            my $len = $subEnd - $subStart + 1;
            #print join(" - ", $offset, $subStart-1, $subEnd)."\n";
            $covStr .= substr($rs->{$acc}->{$offset}, $subStart-1, $len);
        }
        die("$totalLen != ".length($covStr)."\n") if $totalLen != length($covStr);
        
        my @covAry;
        for my $i (0..length($covStr)-1) {
            my $chr = substr($covStr, $i, 1);
            push (@covAry, ord($chr)-33);
        }
        $covHash1->{$acc} = \@covAry;
        
        my $coveredLen = 0;
        for my $cov (@covAry) {
            $coveredLen ++ if $cov >= $cutOffCov;
        }
        $covHash2->{$acc} = $coveredLen/$totalLen;
    }
    return ($covHash1, $covHash2);
}
=usage
my ($covAcc) = getCov(-acc=>["HM005", "HM006"], -chr=>2, -start=>0, -end=>100);
for my $acc (sort(keys(%$covAcc))) {
    my $cov = $covAcc->{$acc};
    my @tmp;
    for my $start (sort(keys(%$cov))) {
        push(@tmp, join("-", $start, $cov->{$start}->[0], $cov->{$start}->[1]));
    }
    print $acc."\t".join(" | ", @tmp), "\n";
}
=cut
sub getCov_old {
    my ($accAry, $chr, $start, $end, $cutOffCov) = rearrange(["acc", "chr", "start", "end", 'cutoffcov'], @_);
    $cutOffCov = $cutOffCov ? $cutOffCov : 1;
    my $chrName = $chr=~/^\d$/ ? "MtChr".$chr : $chr;
    my $totalLen = $end - $start + 1;
    my ($rst1, $rst2) = ({}, {});
    $accAry = ref($accAry) ? $accAry : [$accAry];
    for my $acc (@$accAry) {
        my $covHash;
        my $covLen = 0;
        my $dir = dir($DIR_Mt, "Coverage", $acc);
        my $fIn = dir($dir, $chrName."_stat.txt");
        if(-e $fIn) {
            my $fInH = new IO::File $fIn, "r";
            while(<$fInH>) {
                chomp;
                my ($c, $a, $b) = split("\t", $_);
                next if $b < $start;
                last if $a > $end;
                my $startL = $start<=$a ? $a : $start;
                my $endL   = $end  >=$b ? $b : $end;
                $covHash->{$startL} = [$endL, $c];
                $covLen += $endL - $startL + 1 if $c >= $cutOffCov;
            }
        } else {
            print "no coverage data for $chrName - $acc\n";
            $covHash->{$start} = [$end, 1];
        }
        $rst1->{$acc} = $covHash;
        $rst2->{$acc} = $covLen/$totalLen;
    }
    return ($rst1, $rst2);
}
=usage
my $fIn = file($DIR_Annotation, "Mt3.0_a_BAC_o.gff");
my $fOut = file($DIR_Out, "mt_bac.txt");
chrBand(in=>$fIn, out=>$fOut);
=cut
sub chrBand_old {
    my ($fIn, $fOut) = rearrange(['in', 'out'], @_);
    my $fInH  = new IO::File $fIn, "r";
    my $fOutH = new IO::File $fOut, "w";
    my %phaseDict;
    my %phaseColor = (1=>"grey", 2=>"dgrey", 3=>"vdgrey");
    %phaseColor = (1=>"vlred", 2=>"lred", 3=>"red");
    while(<$fInH>) {
        chomp;
        next if(!$_);
        my @eleAry = split("\t", $_);
        my ($chrName, $source, $type, $start, $stop, $score, $strand, $phase, $desc) =
            ($eleAry[0], $eleAry[1], $eleAry[2], $eleAry[3], $eleAry[4], $eleAry[5], $eleAry[6],
              $eleAry[7], $eleAry[8]);
        my @attrAry1 = split(";", $desc);
        my @attrAry;
        for my $attr (@attrAry1) {
            my @pair = split("=", $attr);
            push (@attrAry, "-".$pair[0], $pair[1]);
        }
        my ($id, $note, $parent) = rearrange(["ID", "Note", "Parent"], @attrAry);
        if($type eq "chromosome") {
            print $fOutH join("\t", "chr", "-", $chrName, $chrName, $start-1, $stop, lc(substr($chrName,2)))."\n";
        } elsif($type eq "BAC") {
            $note =~ /phase\((\d)\)/;
            assert($1);
            $phaseDict{$id} = $1;
        } else {
            assert($type eq "tilingBAC" && $parent);
            $phase = $phaseDict{$parent};
            print $fOutH join("\t", "band", $chrName, $parent, $parent, $start-1, $stop, $phaseColor{$phase})."\n";
        }
    }
}
return 1;
