package VntFilter;
use strict; use Init; use Common; use Localdb; use Readfile; use Parser;
use Path::Class; use DBI; use IO::File; use Vnt; use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Math::Round qw/round/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/applyFilters filter_poss filter_inds filter_snpstat/;
@EXPORT_OK = qw//;
sub filter_inds {
    my ($vnt, $indsN) = @_;
    my ($nInd, $nPos, $inds, $poss, $refs, $snps) = map {$vnt->{$_}} qw/nInd nPos inds poss refs snps/;
    my $vntR = { map {$_ => $vnt->{$_}} qw/nPos poss refs/ };
    my %indH = map {$inds->[$_] => $_} (0..$nInd-1);
    my @indIdxs;
    for my $indN (@$indsN) {
        die "$indN not found\n" unless exists $indH{$indN};
        my $indIdx = $indH{$indN};
        push @indIdxs, $indIdx;
    }
    $vntR->{inds} = [ map {$inds->[$_]} @indIdxs ];
    $vntR->{nInd} = scalar @indIdxs;
    $vntR->{snps} = [ map {$snps->[$_]} @indIdxs ];
    return $vntR;
}
sub filter_poss {
    my ($vnt, $possN) = @_;
    my ($nInd, $nPos, $inds, $poss, $refs, $snps) = map {$vnt->{$_}} qw/nInd nPos inds poss refs snps/;
    my $vntR = { map {$_ => $vnt->{$_}} qw/nInd inds/ };
    my %posH = map {$poss->[$_] => $_} (0..$nPos-1);
    my @posIdxs;
    for my $posN (@$possN) {
        my $posIdx = $posH{$posN};
        die "$posN not found\n".join(" ", @$poss)."\n" unless defined $posIdx;
        push @posIdxs, $posIdx;
    }
    $vntR->{poss} = [ map {$poss->[$_]} @posIdxs ];
    $vntR->{nPos} = scalar @posIdxs;
    $vntR->{refs} = [ map {$refs->[$_]} @posIdxs ];
    my @snpsN;
    for my $snp (@$snps) {
        my @tmp = map {$snp->[$_]} @posIdxs;
        push @snpsN, \@tmp;
    }
    $vntR->{snps} = \@snpsN;
    return $vntR;
}
sub getSnpStat {
    my ($ntRef) = @_;
    my @nts = values %$ntRef;
    my ($flag, $gt, $maf, $acnt);
    my $ntCnt;
    for (@nts) {
        $ntCnt->{$_} ||= 0;
        $ntCnt->{$_} ++;
    }
    my $ntStr = join("", keys %$ntCnt);
    if($ntStr =~ /[^ATCGN]/i) {
        $flag = "ambig";
    } elsif(exists $ntRef->{HM101} && $ntRef->{HM101} eq "N") {
        $flag = "hetero";
    } else {
        $flag = "snp";
    }
    delete $ntCnt->{N} if exists $ntCnt->{N};
    $acnt = scalar(keys %$ntCnt);
    if($acnt == 0) {
        ($gt, $maf) = (0, 0);
    } else { 
        $gt = sum(values %$ntCnt) / @nts;
        $maf = min(values %$ntCnt) / @nts;
    }
    my $h = {flag=>$flag, gt=>$gt, maf=>$maf, acnt=>$acnt};
    return $h;
}
sub filter_snpstat {
    my ($vnt, $co) = @_;
    my ($nInd, $nPos, $inds, $poss, $refs, $snps) = map {$vnt->{$_}} qw/nInd nPos inds poss refs snps/;
    my @poss;
    for my $i (0..$nPos-1) {
        my $pos = $poss->[$i];
        my %nts = map {$inds->[$_] => $snps->[$_]->[$i]} (0..$nInd-1);
        my $stat = getSnpStat(\%nts);
        my ($flag, $gt, $maf, $acnt) = map {$stat->{$_}} qw/flag gt maf acnt/;
        if($flag eq "ambig" || $flag eq "hetero") {
#      printf join("\t", qw/%d %.03f %.03f %d %s/), $pos, $gt, $maf, $acnt, join(" ", @nts)."\n";
        }
        if($flag eq "snp" && $acnt == 2 && $gt >= $co->{gt} && $maf >= $co->{maf}) {
            push @poss, $pos;
        }
    }
    return filter_poss($vnt, \@poss);
}
sub applyFilters {
    my ($vnt, $inds, $poss, $co) = rearrange(['vnt', 'inds', 'poss', 'co'], @_);
    printf "input:          %3d inds, %6d poss\n", map {$vnt->{$_}} qw/nInd nPos/;
    $vnt = filter_inds($vnt, $inds) if $inds && @$inds > 0;
    printf "ind_filter:     %3d inds, %6d poss\n", map {$vnt->{$_}} qw/nInd nPos/;
    $vnt = filter_poss($vnt, $poss) if $poss && @$poss > 0;
    printf "pos_filter:     %3d inds, %6d poss\n", map {$vnt->{$_}} qw/nInd nPos/;
    $vnt = filter_snpstat($vnt, $co) if $co;
    printf "snpstat_filter: %3d inds, %6d poss\n", map {$vnt->{$_}} qw/nInd nPos/;
    return $vnt;
}


1;
__END__
