#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal.calib.pl - Fills the stat fields (match, misMatch, baseN)

=head1 SYNOPSIS
  
  gal.calib.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gal)
    -o (--out)    output file (Gal)
    -q (--qry)    query fasta
    -t (--tgt)    target fasta

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Time::HiRes qw/gettimeofday tv_interval/;
use Location;
use Seq;
use Gal;

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

my $t0 = [gettimeofday];
print $fho join("\t", @HEAD_GAL)."\n";

my $cnt = 1;
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps == 21;
  my ($cid, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $lev, $ali, $mat, $mis, $qN, $tN, $ident, $score, $tlS, $qlS) = @$ps;
  my ($rqLoc, $rtLoc) = (locStr2Ary($qlS), locStr2Ary($tlS));
  @$rqLoc == @$rtLoc || die "unequal pieces\n";
  my ($rqb, $rqe) = ($rqLoc->[0]->[0], $rqLoc->[-1]->[1]);
  my ($rtb, $rte) = ($rtLoc->[0]->[0], $rtLoc->[-1]->[1]);
  $qEnd - $qBeg == $rqe - $rqb ||
    die "qLen err: $cid $qId($qBeg-$qEnd): $rqb-$rqe\n"; 
  $tEnd - $tBeg == $rte - $rtb ||
    die "tLen err: $cid $tId($tBeg-$tEnd): $rtb-$rte\n";
  my ($qlen, $tlen) = (locAryLen($rqLoc), locAryLen($rtLoc));
  $qlen == $tlen || die "unequal loclen: qry[$qlen] - tgt[$tlen]\n";
  
  my $nBlock = @$rqLoc;

#    my $seqT = seqRet([[$tBeg, $tEnd]], $tId, $tSrd, $ft);
#    my $seqQ = seqRet([[$qBeg, $qEnd]], $qId, $qSrd, $fq);
#    my $tSeq = getSubSeq($seqT, $rtLoc);
#    my $qSeq = getSubSeq($seqQ, $rqLoc);
  my $tLoc = $tSrd eq "-" ? 
    [ map {[$tEnd-$_->[1]+1, $tEnd-$_->[0]+1]} @$rtLoc ]
    : [ map {[$tBeg+$_->[0]-1, $tBeg+$_->[1]-1]} @$rtLoc ]; 
  my $qLoc = $qSrd eq "-" ? 
    [ map {[$qEnd-$_->[1]+1, $qEnd-$_->[0]+1]} @$rqLoc ]
    : [ map {[$qBeg+$_->[0]-1, $qBeg+$_->[1]-1]} @$rqLoc ];
  
  my $tSeq = seqRet($tLoc, $tId, $tSrd, $ft);
  my $qSeq = seqRet($qLoc, $qId, $qSrd, $fq);
  if(length($tSeq) != $tlen) {
    print "tseq error\n$tId:$tBeg-$tEnd $tlS\n$tlen\n".length($tSeq)."\n";
    exit;
  }
  if(length($qSeq) != $qlen) {
    print "tseq error\n$qId:$qBeg-$qEnd $qlS\n$qlen\n".length($qSeq)."\n";
    exit;
  }
  ($mat, $mis, $qN, $tN) = seqCompare($tSeq, $qSeq);
  $ali = $mat + $mis + $qN + $tN;
  $ident = ($mat+$mis==0) ? 0 : sprintf "%.03f", $mat/($mat+$mis);
  @$ps[12..17] = ($ali, $mat, $mis, $qN, $tN, $ident);
  print $fho join("\t", @$ps)."\n";

  printf "%5d: %.01f min\n", $cnt++, 
    tv_interval($t0, [gettimeofday]) / 60 if $cnt % 1000 == 0;
}
close $fhi;
close $fho;


__END__
