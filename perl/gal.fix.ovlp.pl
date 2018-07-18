#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal.fix.ovlp.pl - fix ovlp-blocks in a Gal file

=head1 SYNOPSIS
  
  gal.fix.ovlp.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file

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
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

print $fho join("\t", @HEAD_GAL)."\n";

my $cnt = 0;
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps == 21;
  my ($cid, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $lev, $ali, $mat, $mis, $qN, $tN, $ident, $score, $tlS, $qlS) = @$ps;
  $tSrd eq "+" || die "$cid: tSrd -\n";
  
  my ($qLoc, $tLoc) = (locStr2Ary($qlS), locStr2Ary($tlS));
  my @qlens = map {$_->[1]-$_->[0]+1} @$qLoc;
  my $ref = tiling($qLoc, \@qlens, 2);
  my (@rqloc, @rtloc);
  for (@$ref) {
    my ($rqb, $rqe, $idx) = @$_;

    my ($qb, $qe) = @{$qLoc->[$idx]};
    my ($tb, $te) = @{$tLoc->[$idx]};
    my $rtb = $rqb - $qb + $tb;
    my $rte = $rqe - $qb + $tb;

    if(@rqloc == 0 || $rtb > $rtloc[-1]->[1]) { 
        push @rqloc, [$rqb, $rqe];
        push @rtloc, [$rtb, $rte];
    }
  }
  my ($nqLocS, $ntLocS) = (locAry2Str(\@rqloc), locAry2Str(\@rtloc));
  if($nqLocS ne $qlS) {
    @$ps[19,20] = ($ntLocS, $nqLocS);
    $cnt ++;
  }
  print $fho join("\t", @$ps)."\n";
}
print STDERR "$cnt rows fixed\n";
close $fhi;
close $fho;


__END__
