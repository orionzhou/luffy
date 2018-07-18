#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal2gax.pl - convert a Gal (wide) file to Gax (long) file

=head1 SYNOPSIS
  
  gal2gax.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gal)
    -o (--out)    output file (Gax: tid tb te tsrd qid qb qe qsrd cid lev)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use Gal;

my ($fi, $fo) = ('') x 2;
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

while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps == 21;
  my ($cid, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $lev, $ali, $mat, $mis, $qN, $tN, $ident, $score, $tlS, $qlS) = @$ps;
  my ($rqloc, $rtloc) = (locStr2Ary($qlS), locStr2Ary($tlS));
  for my $i (0..@$rqloc-1) {
    my ($rtb, $rte) = @{$rtloc->[$i]};
    my ($rqb, $rqe) = @{$rqloc->[$i]};
    my ($tb, $te) = $tSrd eq "-" ? ($tEnd-$rte+1, $tEnd-$rtb+1)
      : ($tBeg+$rtb-1, $tBeg+$rte-1);
    my ($qb, $qe) = $qSrd eq "-" ? ($qEnd-$rqe+1, $qEnd-$rqb+1)
      : ($qBeg+$rqb-1, $qBeg+$rqe-1);
    print $fho join("\t", $tId, $tb, $te, $tSrd, 
      $qId, $qb, $qe, $qSrd, $cid, $lev)."\n";
  }
}
close $fhi;
close $fho;


__END__
