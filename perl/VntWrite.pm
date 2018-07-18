package VntWrite;
use strict; use Init; use Common; use Localdb; use Readfile;
use Path::Class; use DBI; use IO::File; use Switch; use Data::Dumper; 
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/vntWrite/;
@EXPORT_OK = qw//;
sub vntWrite {
    my ($format, $pre, $vnt) = rearrange(["format", "pre", "vnt"], @_);
    if($format eq "simplesnp") {
        vnt_write_simplesnp($vnt, $pre);
    } elsif($format eq "haploview") {
        vnt_write_haploview($vnt, $pre);
    } elsif($format eq "phase") { 
        vnt_write_phase($vnt, $pre);
    } elsif($format eq "fastphase") { 
        vnt_write_fastphase($vnt, $pre);
    } elsif($format eq "rsq") { 
        vnt_write_rsq($vnt, $pre);
    } elsif($format eq "ldhat") { 
        vnt_write_ldhat($vnt, $pre);
    } elsif($format eq "ldhot") {
        vnt_write_ldhot($vnt, $pre);
    } else {
        die("format[$format] not supported\n") unless $format eq "table";
        vnt_write_table($vnt, $pre);
    }
}
sub vnt_write_simplesnp {
    my ($vnt, $pre) = @_;
    my ($nInd, $nPos, $inds, $poss, $refs, $snps) = map {$vnt->{$_}} qw/nInd nPos inds poss refs snps/;
    my $fh = new IO::File "$pre.txt", "w";
    print $fh join(" ", $nInd, $nPos)."\n";
    print $fh join(" ", @$poss)."\n";
    print $fh join(" ", @$refs)."\n";
    for my $i (0..$nInd-1) {
        my ($ind, $snp) = ($inds->[$i], $snps->[$i]);
        my $snpStr = join(" ", @$snp);
        print $fh join(" ", $ind, $snpStr)."\n";;
    }
}
sub vnt_write_haploview {
    my ($vnt, $pre) = @_;
    my ($nInd, $nPos, $inds, $poss, $refs, $snps) = map {$vnt->{$_}} qw/nInd nPos inds poss refs snps/;
    my $fh1 = new IO::File "$pre.txt", "w";
    my $fh2 = new IO::File "$pre\_loc.txt", "w";
    for my $i (0..$nInd-1) {
        my ($ind, $snp) = ($inds->[$i], $snps->[$i]);
        my $snpStr = join(" ", @$snp);
        $snpStr =~ tr/NACGT/01234/;
        print $fh1 join("\t", "Mt", $ind, $snpStr)."\n";;
    }
    for my $j (0..($nPos-1)) {
        my $pos = $poss->[$j];
        print $fh2 join("\t", sprintf("L%0".getDigits($nPos)."d", $j+1), $pos)."\n"; 
    }
}
sub vnt_write_phase {
    my ($vnt, $pre) = @_;
    my ($nInd, $nPos, $inds, $poss, $refs, $snps) = map {$vnt->{$_}} qw/nInd nPos inds poss refs snps/;
    my $fh = new IO::File "$pre.txt", "w";
    print $fh join("\n", $nInd, $nPos, "S"x$nPos, join(" ", @$poss))."\n";
    for my $i (0..$nInd-1) {
        my ($ind, $snp) = ($inds->[$i], $snps->[$i]);
        print $fh "$ind\n";
        my $snpStr = join("", @$snp);
        $snpStr =~ s/N/\?/g;
        print $fh join("\n", ($snpStr) x 2)."\n";
    }
}
sub vnt_write_fastphase {
    my ($vnt, $pre) = @_;
    my ($nInd, $nPos, $inds, $poss, $refs, $snps) = map {$vnt->{$_}} qw/nInd nPos inds poss refs snps/;
    my $fh = new IO::File "$pre.txt", "w";
    print $fh join("\n", $nInd, $nPos, join(" ", "P", @$poss), "S"x$nPos)."\n";
    for my $i (0..$nInd-1) {
        my ($ind, $snp) = ($inds->[$i], $snps->[$i]);
        print $fh "$ind\n";
        my $snpStr = join("", @$snp);
        $snpStr =~ s/N/\?/g;
        print $fh join("\n", ($snpStr) x 2)."\n";
    }
}
sub vnt_write_rsq {
    my ($vnt, $pre) = @_;
    my ($nInd, $nPos, $inds, $poss, $refs, $snps) = map {$vnt->{$_}} qw/nInd nPos inds poss refs snps/;
    my $fh = new IO::File "$pre.txt", "w";
    print $fh join("\n", $nInd, $nPos, join(" ", "P", @$poss), "S"x$nPos)."\n";
    for my $i (0..$nInd-1) {
        my ($ind, $snp) = ($inds->[$i], $snps->[$i]);
        my $snpStr = join("", @$snp);
        print $fh ">$ind\n";
        print $fh "$snpStr\n";
    }
}
sub vnt_write_ldhat {
    my ($vnt, $pre) = @_;
    my ($nInd, $nPos, $inds, $poss, $refs, $snps) = map {$vnt->{$_}} qw/nInd nPos inds poss refs snps/;
    my $fh1 = new IO::File "$pre.txt", "w";
    my $fh2 = new IO::File "$pre\_loc.txt", "w";
    my ($posMin, $posMax) = ($poss->[0], $poss->[-1]);
    print $fh2 sprintf("%d\t%d\tL", $nPos, ceil(($posMax-$posMin+1)/1000))."\n"; #LDhat_loc
    for my $i (0..($nPos-1)) {
        my $pos = $poss->[$i];
        print $fh2 sprintf("%.3f", ($pos-$posMin+1)/1000)."\n";
    }
    print $fh1 sprintf("%d %d 1", $nInd, $nPos)."\n"; 
    for my $i (0..$nInd-1) {
        my ($ind, $snp) = ($inds->[$i], $snps->[$i]);
        my $snpStr = join("", @$snp);
        $snpStr =~ s/N/\?/g;
        print $fh1 ">$ind\n";
        print $fh1 "$snpStr\n";
    }
}
sub vnt_write_ldhot {
    my ($vnt, $pre) = @_;
    my ($nInd, $nPos, $inds, $poss, $refs, $snps) = map {$vnt->{$_}} qw/nInd nPos inds poss refs snps/;
    my $fh = new IO::File "$pre.txt", "w";
    my ($posMin, $posMax) = ($poss->[0], $poss->[-1]);
    print $fh "Distinct = $nInd\n"
                ."Genes = $nInd\n"
                ."Loci = $nPos\n"
                ."K = -4\n"
                ."Positions of loci:\n"
                .join("\t", @$poss)."\n"
                ."Haplotypes\n";
    for my $i (0..$nInd-1) {
        my ($ind, $snp) = ($inds->[$i], $snps->[$i]);
        my $snpStr = join("", @$snp);
        $snpStr =~ s/N/\?/g;
        print $fh "\t$snpStr\t1\n";
    }
    print $fh "#\n";
}
sub vnt_write_table {
    my ($vnt, $pre) = @_;
    my ($nInd, $nPos, $inds, $poss, $refs, $snps) = map {$vnt->{$_}} qw/nInd nPos inds poss refs snps/;
    my $fh = new IO::File "$pre.txt", "w";
    print $fh join("\t", "position", @$inds)."\n";
    for my $j (0..$nPos-1) {
        my @alleles = map {$snps->[$_]->[$j]} (0..$nInd-1);
        print $fh join("\t", $poss->[$j], prettyStr(join("", @alleles)))."\n";
    }
}



1;
__END__
