#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal2idm.pl - Call Ins/Del/Mnp from a Gal file

=head1 SYNOPSIS
  
  gal2idm.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gal)
    -o (--out)    output tabular file (tid tb te qid qb qe cid lev)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use Common;
use Seq;
use File::Path qw/make_path remove_tree/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fi, $fo) = ('') x 2;
my ($len) = (50); 
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "len|l=i"  => \$len,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

open ($fho, ">$fo") || die "cannot write $fo\n";

while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps == 21;
  my ($cid, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $lev, $ali, $mat, $mis, $qN, $tN, $ident, $score, $tlS, $qlS) = @$ps;
  my ($rqLoc, $rtLoc) = (locStr2Ary($qlS), locStr2Ary($tlS));
  
  $rqLoc = [ sort {$a->[0] <=> $b->[0]} @$rqLoc ];
  $rtLoc = [ sort {$a->[0] <=> $b->[0]} @$rtLoc ];
  @$rqLoc == @$rtLoc || die "unequal pieces\n";
  my $nBlock = @$rqLoc;
  my @lens = map {$_->[1] - $_->[0] + 1} @$rqLoc;
 
  $nBlock > 1 || next;
  for my $i (1..$nBlock-1) {
    my ($rqb, $rqe) = ($rqLoc->[$i-1]->[1], $rqLoc->[$i]->[0]);
    my ($rtb, $rte) = ($rtLoc->[$i-1]->[1], $rtLoc->[$i]->[0]);
    my $qins = $rqe - $rqb - 1;
    my $tins = $rte - $rtb - 1;
    my ($qb, $qe) = $qSrd eq "-" ? ($qEnd - $rqe + 1, $qEnd - $rqb + 1) : 
      ($qBeg + $rqb - 1, $qBeg + $rqe - 1);
    my ($tb, $te) = $tSrd eq "-" ? ($tEnd - $rte + 1, $tEnd - $rtb + 1) : 
      ($tBeg + $rtb - 1, $tBeg + $rte - 1);
    print $fho join("\t", $tId, $tb, $te, $tSrd, 
      $qId, $qb, $qe, $qSrd, $cid, $lev)."\n";
  }
}
close $fhi;
close $fho;


__END__
