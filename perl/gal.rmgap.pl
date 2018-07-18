#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal.rmgap.pl - remove assembly gaps within alignment blocks

=head1 SYNOPSIS
  
  gal.rmgap.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -q (--qry)    qry (BED) file with gap locations
    -t (--tgt)    tgt (BED) file with gap locations 

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Data::Dumper;
use Location;
use Gal;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fi, $fo) = ('') x 2;
my ($fq, $ft) = ('') x 2; 
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "qry|q=s"  => \$fq,
  "tgt|t=s"  => \$ft,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fq || !$ft;

if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $hq = get_gap_hash($fq);
my $ht = get_gap_hash($ft);

print $fho join("\t", @HEAD_GAL)."\n";

while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps == 21;
  my ($cid, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $lev, $ali, $mat, $mis, $qN, $tN, $ident, $score, $tlS, $qlS) = @$ps;
  my ($rqLoc, $rtLoc) = (locStr2Ary($qlS), locStr2Ary($tlS));
  my ($orqLoc, $ortLoc) = ($rqLoc, $rtLoc);
  
#  next unless $id == 32314;
  my $qgLoc = $hq->{$qId};
  my $rqgLoc = $qSrd eq "+" ? 
    [ map {[$_->[0]-$qBeg+1, $_->[1]-$qBeg+1]} @$qgLoc ] :
    [ map {[$qEnd-$_->[1]+1, $qEnd-$_->[0]+1]} @$qgLoc ];
  my ($qoLoc, $qoLen) = posOvlp($rqgLoc, $rqLoc);
  if($qoLen > 0) {
    my ($nrqLoc) = posDiff($rqLoc, $qoLoc);
    my $toLoc = coordMap($rqLoc, $rtLoc, "+", "+", $qoLoc);
    my ($nrtLoc) = posDiff($rtLoc, $toLoc);
    $rqLoc = $nrqLoc;
    $rtLoc = $nrtLoc;
  }

  my $tgLoc = $ht->{$tId};
  my $rtgLoc = $tSrd eq "+" ? 
    [ map {[$_->[0]-$tBeg+1, $_->[1]-$tBeg+1]} @$tgLoc ] :
    [ map {[$tEnd-$_->[1]+1, $tEnd-$_->[0]+1]} @$tgLoc ];
  my ($toLoc, $toLen) = posOvlp($rtgLoc, $rtLoc);
  if($toLen > 0) {
    my ($nrtLoc) = posDiff($rtLoc, $toLoc);
    my $qoLoc = coordMap($rtLoc, $rqLoc, "+", "+", $toLoc);
    my ($nrqLoc) = posDiff($rqLoc, $qoLoc);
    $rqLoc = $nrqLoc;
    $rtLoc = $nrtLoc;
  }
  
  my $osb = $rqLoc->[0]->[0] - 1;
  my $ose = $orqLoc->[-1]->[1] - $rqLoc->[-1]->[1];
  if($qSrd eq "+") {
    $qBeg += $osb if $osb != 0;
    $qEnd -= $ose if $ose != 0;
  } else {
    $qEnd -= $osb if $osb != 0;
    $qBeg += $ose if $ose != 0;
  }
  if($tSrd eq "+") {
    $tBeg += $osb if $osb != 0;
    $tEnd -= $ose if $ose != 0;
  } else {
    $tEnd -= $osb if $osb != 0;
    $tBeg += $ose if $ose != 0;
  }
  @$ps[2,3,7,8] = ($tBeg, $tEnd, $qBeg, $qEnd);
  if($osb != 0) {
    $rqLoc = [ map {[$_->[0] - $osb, $_->[1] - $osb]} @$rqLoc ];
    $rtLoc = [ map {[$_->[0] - $osb, $_->[1] - $osb]} @$rtLoc ];
  }
#  if($id == 32314) {
#    print join("\t", $qoLen, $toLen, $osb, $ose)."\n";
#    print join("\t", locAryLen($rqLoc), locAryLen($rtLoc))."\n";
#    die;
#  }
  $ps->[19] = locAry2Str($rtLoc);
  $ps->[20] = locAry2Str($rqLoc);
  print $fho join("\t", @$ps)."\n";
}
close $fhi;
close $fho;

sub get_gap_hash {
  my ($fi) = @_;
  my $h = {};
  if(-e $fi && -s $fi) {
    my $ti = readTable(-in=>$fi, -header=>0);
    for my $i (0..$ti->lastRow) {
      my ($chr, $beg, $end) = $ti->row($i);
      $h->{$chr} ||= [];
      push @{$h->{$chr}}, [$beg + 1, $end];
    }
  } elsif(-e $fi && ! -s $fi) {  ## do nothing
  } else {
    die "$fi is not there\n";
  }
  return $h;
}
sub coordMap {
  my ($qloc, $tloc, $qsrd, $tsrd, $iloc) = @_;
  $qloc = [ sort {$a->[0] <=> $b->[0]} @$qloc ];
  $qloc = [ reverse @$qloc ] if $qsrd =~ /^\-1?$/;
  $tloc = [ sort {$a->[0] <=> $b->[0]} @$tloc ];
  $tloc = [ reverse @$tloc ] if $tsrd =~ /^\-1?$/;
  my $opp = is_revsrd($qsrd, $tsrd);

  @$qloc == @$tloc || 
    die "unequal loc pieces\n".Dumper($qloc).Dumper($tloc);
  my @oloc;
  for my $i (0..@$iloc-1) {
    my ($ib, $ie) = @{$iloc->[$i]};
    my $idx = first_index {$_->[0] <= $ib && $ie <= $_->[1]} @$qloc;
    $idx > -1 || die "[$ib-$ie] not in locI\n".Dumper($qloc);
    my ($qb, $qe) = @{$qloc->[$idx]};
    my ($tb, $te) = @{$tloc->[$idx]};
    my $ob = $opp ? $te-($ib-$qb) : $tb+($ib-$qb);
    my $oe = $opp ? $te-($ie-$qb) : $tb+($ie-$qb);
    push @oloc, [$ob, $oe];
  }
  return \@oloc;
}

__END__
