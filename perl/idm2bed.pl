#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  idm2bed.pl - convert Idm file to Bed file

=head1 SYNOPSIS
  
  idm2bed.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input (IDM) file
      -out    output (BED) file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
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
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

while( <$fhi> ) {
  chomp;
  next if /(^\#)|(^\s*$)/;
  my ($tid, $tb, $te, $ts, $qid, $qb, $qe, $qs, $id, $lev) = split "\t";
  my $tlen = $te - $tb + 1;
  my $qlen = $qe - $qb + 1;
  my $name = "$tlen/$qlen";
  my $score = exists $h_score->{$lev} ? $h_score->{$lev} : 300;
  print $fho join("\t", $tid, $tb - 1, $te - 1, $name, $score)."\n";
}
close $fhi;
close $fho;


__END__
