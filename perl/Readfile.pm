package Readfile;
use strict;
use Init;
use Common;
use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/readLocStr getIds getAccInDir
    getCov getRegionCov getVnt
    getSeqDesc getSeqInfo/;
@EXPORT_OK = qw//;
sub getIds {
    my ($fi, $idCol) = @_;
    die "$fi is not there\n" unless -s $fi;
    $idCol ||= "id";

    my $t = readTable(-in=>$fi, -header=>1);
    my @ids = grep /\S+/, $t->col($idCol);
    @ids = uniq(@ids);
    for (@ids) {
        $_ =~ s/(^IMGAG?\|)|(\|PAC.*$)//i;
    }
    printf "\t%5d ids obtained from $fi\n", scalar(@ids);
    return \@ids;
}
sub getCov {
    my ($loc, $covOpt, $refDb) = rearrange(["loc", "covopt", "refdb"], @_);
    $covOpt ||= 1;
    my ($seqid, $beg, $end) = ($loc->seq_id, $loc->start, $loc->end);
    my $len = $end - $beg + 1;
    die "not uniform seqid\n".Dumper($loc) unless $seqid;

    my $fi;
    if($refDb eq "mt_35") {
        $fi = "$DIR_Misc3/hapmap/30_vnt/21_coverage";
    } else {
        die "unknonw refdb: $refDb\n" unless $refDb eq "mt_30";
        $fi = "$DIR_Misc3/hapmap_mt30/30_vnt/21_coverage";
    }
    if($covOpt == 1) {
        $fi = "$fi/cov.tbl.gz";
    } else {
        die "unknown opt: $covOpt\n" unless $covOpt == 2;
        $fi = "$fi/cov_uniq.tbl.gz";
    }
    die "$fi is not there\n" unless -s $fi;

    my @ids_all = $refDb eq "mt_30" ? 
        map {sprintf "HM%03d", $_} (1..30,101) :
        map {sprintf "HM%03d", $_} (1..60,101,102,318..339);
    my $idx_offset = 3;
    my $idH = { map {$ids_all[$_] => $_+$idx_offset} (0..$#ids_all) };

    my $region = sprintf "%s:%d-%d", $loc->seq_id, $loc->start, $loc->end;
    my $cmd = "tabix $fi $region";
    open(my $fh, $cmd." |") || die "Failed: $! in \n$cmd\n";
    return sub {
        return undef if eof($fh);
        my $line = readline($fh);
        chomp $line;
        my @ps = split("\t", $line);
        my ($chr, $pos) = @ps[0..1];
        my $posR = $pos - $beg + 1;
        my $ref_cov = { map {$_ => $ps[$idH->{$_}] ? $ps[$idH->{$_}] : 0} @ids_all };
        my $ref = {chr=>$chr, pos=>$pos, posR=>$posR, cov=>$ref_cov};
        return $ref;
    }
}
sub getVnt {
    my ($loc, $refDb) = rearrange(["loc", "refdb"], @_);
    my ($seqid, $beg, $end) = ($loc->seq_id, $loc->start, $loc->end);
    my $len = $end - $beg + 1;
    die "not uniform seqid\n".Dumper($loc) unless $seqid;

    my $fi;
    if($refDb eq "mt_35") {
        $fi = "$DIR_Misc3/hapmap/30_vnt/11_tbl/01.tbl.gz";
    } else {
        die "unknonw refdb: $refDb\n" unless $refDb eq "mt_30";
        $fi = "$DIR_Misc3/hapmap_mt30/30_vnt/11_tbl/01.tbl.gz";
    }
    die "$fi is not there\n" unless -s $fi;

    my @ids_all = $refDb eq "mt_30" ? 
        map {sprintf "HM%03d", $_} (1..30,101) :
        map {sprintf "HM%03d", $_} (1..60,101,102,318..339);
    my $idx_offset = 4;
    my $idH = { map {$ids_all[$_] => $_+$idx_offset} (0..$#ids_all) };

    my $region = sprintf "%s:%d-%d", $seqid, $beg, $end;
    my $cmd = "tabix $fi $region";
    open(my $fh, $cmd." |") || die "Failed: $! in \n$cmd\n";
    return sub {
        return undef if eof($fh);
        my $line = readline($fh);
        chomp $line;
        my @ps = split("\t", $line);
        my ($chr, $pos, $refa) = @ps[0..1, 3];
        my $posR = $pos - $beg + 1;
        my $vnt = { map {$_ => $ps[$idH->{$_}]} @ids_all };
        my $ref = {chr=>$chr, pos=>$pos, refa=>$refa, posR=>$posR, vnt=>$vnt};
        return $ref;
    }
}
sub getRegionCov {
    my ($loc, $covOpt, $accs, $refDb) = @_;
    my $it = getCov(-loc=>$loc, -covopt=>$covOpt, -refdb=>$refDb);
    my $covh = { map {$_=>[]} @$accs };
    while( my $ref = $it->() ) {
        my ($chr, $pos, $posR) = map {$ref->{$_}} qw/chr pos posR/;
        for my $acc (@$accs) {
            push @{$covh->{$acc}}, $ref->{cov}->{$acc};
        }
    }
    return $covh;
}
sub readLocStr {
    my ($fIn, $refDb) = @_;
    unless($refDb) {
        $refDb = 'mt_30';
        print "using $refDb as refDb\n";
    }
    if(!$fIn || ! -s $fIn) {
        $fIn = file($DIR_In, 'mt.txt');
        print "locFile using $fIn\n";
    }
    my $locH = {};
    die "$fIn is not there\n" unless $fIn;
    my $fInH = new IO::File $fIn, "r";
    while( <$fInH> ) {
        chomp;
        last if $_ =~ /^--/;
        next unless $_;
        my ($seqid, $start, $end, $strand) = splitLocStr($_);
        $end = min($end, getSeqLen($seqid, $refDb));
        $locH->{$seqid} = [] unless exists $locH->{$seqid};
        push @{$locH->{$seqid}}, [$start, $end];
    }
    return $locH;
}
sub getAccInDir {
    my ($dir) = @_;
    opendir(my $in, $dir) || die "Can't open $dir\n";
    my @accAry;
    for my $name ( readdir($in) ) {
        if($name =~ /^(HM\d+)(.*)$/) {
            my $acc = $1;
            my $idx = first_index {$_ eq $acc} @accAry;
            push @accAry, $acc if $idx == -1;
        }
    }
    return \@accAry;
}




1;
__END__
