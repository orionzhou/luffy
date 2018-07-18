package BioFeature;
use strict;
use Bio::SeqFeature::Generic;
use Common;
use Data::Dumper;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/gtb2Feature loc2Feature
    gtb2Features/;
@EXPORT_OK = qw//;

sub loc2Feature {
    my ($locC, $str, $id, $chr, $locE, $loc5, $loc3) = rearrange([qw/locC strand id chr locE loc5 loc3/], @_);
    $str = (defined($str) && $str =~ /^\-1?$/) ? -1 : 1;
    $locE = [] unless defined($locE);
    $loc5 = [] unless defined($loc5);
    $loc3 = [] unless defined($loc3);
    my $beg = min( map {$_->[0]} (@$locC, @$locE, @$loc5, @$loc3) );
    my $end = max( map {$_->[1]} (@$locC, @$locE, @$loc5, @$loc3) );
    my $fe = Bio::SeqFeature::Generic->new(-display_name=>$id, -start=>$beg, -end=>$end, -strand=>$str, -primary_tag=>'mRNA');
    for (@$locC) {
        my ($beg, $end) = @$_;
        $fe->add_SeqFeature(Bio::SeqFeature::Generic->new(-start=>$beg, -end=>$end, -strand=>$str, -primary_tag=>"CDS"));
    }
    for (@$locE) {
        my ($beg, $end) = @$_;
        $fe->add_SeqFeature(Bio::SeqFeature::Generic->new(-start=>$beg, -end=>$end, -strand=>$str, -primary_tag=>"exon"));
    }
    for (@$loc5) {
        my ($beg, $end) = @$_;
        $fe->add_SeqFeature(Bio::SeqFeature::Generic->new(-start=>$beg, -end=>$end, -strand=>$str, -primary_tag=>"five_prime_UTR"));
    }
    for (@$loc3) {
        my ($beg, $end) = @$_;
        $fe->add_SeqFeature(Bio::SeqFeature::Generic->new(-start=>$beg, -end=>$end, -strand=>$str, -primary_tag=>"three_prime_UTR"));
    }
    return $fe;
}
sub gtb2Feature {
    my ($rowRef) = @_;
    my ($id, $str, $locCS, $locES, $loc5S, $loc3S, $cat1) = @$rowRef[0,5,8,6,9,10,14];
    my ($locC, $locE, $loc5, $loc3) = map {locStr2Ary($_)} ($locCS, $locES, $loc5S, $loc3S);
    $id = "$id [$cat1]" if $cat1 ne "gene";
    my $fe = loc2Feature(-id=>$id, -strand=>$str, -locC=>$locC, -locE=>$locE);
    return $fe;
}
sub gtb2Features {
    my ($tg, $chr, $beg, $end) = @_;
    my @fes;
    my $tg_r = $tg->match_pattern("\$_->[2] eq '$chr' &&
      ( (\$_->[3] > $beg && \$_->[3] < $end) || (\$_->[4] > $beg && \$_->[4] < $end) )");
    for my $i (0..$tg_r->nofRow-1) {
        push @fes, gtb2Feature($tg_r->rowRef($i));
    }
    return \@fes;
}



1;
__END__

