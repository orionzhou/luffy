#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal2chain.pl - convert GAL file to Chain format

=head1 SYNOPSIS
  
  gal2chain.pl [-help] [-in input-file] [-out output-file]

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
use List::Util qw/min max sum/;

my ($fi, $fo) = ('', '');
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

while( <$fhi> ) {
  chomp;
  /(^id)|(^\#)|(^\s*$)/ && next;
  my $ps = [ split "\t" ];
  next unless @$ps == 21;
  my ($cid, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $lev, $ali, $mat, $mis, $qN, $tN, $ident, $score, $tlS, $qlS) = @$ps;
  $tSrd eq "+" || die "$cid: tSrd -\n";

  my ($qloc, $tloc) = (locStr2Ary($qlS), locStr2Ary($tlS));
  @$qloc == @$tloc || die "unequal pieces\n";
  my $nBlock = @$qloc;
  
  my ($ctb, $cte) = ($tBeg - 1, $tEnd);
  my ($cqb, $cqe) = $qSrd eq "+" ? ($qBeg - 1, $qEnd) : 
    ($qSize - $qEnd + 1 - 1, $qSize - $qBeg + 1);
  print $fho join(" ", "chain", $score, $tId, $tSize, $tSrd, $ctb, $cte,
    $qId, $qSize, $qSrd, $cqb, $cqe, $cid)."\n";
  if($nBlock > 1) {
    for my $i (0..$nBlock-2) {
      my $len = $tloc->[$i]->[1] - $tloc->[$i]->[0] + 1;
      my $dt = $tloc->[$i+1]->[0] - $tloc->[$i]->[1] - 1;
      my $dq = $qloc->[$i+1]->[0] - $qloc->[$i]->[1] - 1;
      print $fho join("\t", $len, $dt, $dq)."\n";
    }
  }
  my $tlocl = $tloc->[$nBlock-1];
  my $lenl = $tlocl->[1] - $tlocl->[0] + 1;
  print $fho $lenl."\n\n";
}
close $fhi;
close $fho;

__END__
