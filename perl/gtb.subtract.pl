#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
    gtb.subtract.pl - subtract a qry-Gtb from tgt-Gtb

=head1 SYNOPSIS
  
  gtb_compare.pl [-help] [-qry query-Gtb] [-tgt target-Gtb] [-out subtracted-Gtb]

  Options:
    -h (--help)  brief help message
    -q (--qry)   query Gtb file
    -t (--tgt)   target Gtb file
    -o (--out)   output (subtracted) Gtb file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Location;
use Gtb;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fhq, $fht, $fho);
my ($fq, $ft, $fo) = ('') x 4;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "qry|q=s"  => \$fq,
  "tgt|t=s"  => \$ft,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fq || !$ft || !$fo;

open ($fhq, "<$fq") || die "cannot read $fq\n";
open ($fht, "<$ft") || die "cannot read $ft\n";
open ($fho, ">$fo") || die "cannot write $fo\n";

my $tq = readTable(-inh => $fhq, -header => 1);
my $tt = readTable(-inh => $fht, -header => 1);
close $fhq;
close $fht;

my $hq;
for my $i (0..$tq->nofRow-1) {
  my ($chr, $beg, $end, $srd, $locC)  = map {$tq->elm($i, $_)} qw/chr beg end srd locC/;
  my $rloc = [sort {$a->[0] <=> $b->[0]} @{locStr2Ary($locC)}];
  my $loc = $srd eq "-" ? [map {[$end-$_->[1]+1, $end-$_->[0]+1]} @$rloc] : 
    [map {[$beg+$_->[0]-1, $beg+$_->[1]-1]} @$rloc];
  my $str = join("|", $chr, locAry2Str($loc));
  $hq->{$str} = $i;
}

my $ht;
my @idxs_rm;
for my $i (0..$tt->nofRow-1) {
  my ($chr, $beg, $end, $srd, $locC)  = map {$tt->elm($i, $_)} qw/chr beg end srd locC/;
  my $rloc = [sort {$a->[0] <=> $b->[0]} @{locStr2Ary($locC)}];
  my $loc = $srd eq "-" ? [map {[$end-$_->[1]+1, $end-$_->[0]+1]} @$rloc] : 
    [map {[$beg+$_->[0]-1, $beg+$_->[1]-1]} @$rloc];
  my $str = join("|", $chr, locAry2Str($loc));
  if(exists($ht->{$str})) {
    push @idxs_rm, $ht->{$str};
  }
  $ht->{$str} = $i;
  if(exists($hq->{$str})) {
    push @idxs_rm, $i;
  }
}

printf STDOUT "%5d | %5d substracted\n", scalar(@idxs_rm), $tt->nofRow;

$tt->delRows(\@idxs_rm);
print $fho $tt->tsv(1);
close $fho;


__END__
