#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.phytozome.pl - process Phytozome GTB files

=head1 SYNOPSIS
  
  gtb.phytozome.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (*.gtb)
    -o (--out)    output file (*.gtb)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Gtb;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my ($fhi, $fho);
if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot read $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $hs;
my $cnt = 0;
print $fho join("\t", @HEAD_GTB)."\n";
while( <$fhi> ) {
  chomp;
  /(^id)|(^\#)|(^\s*$)/i && next;
  my $ps = [ split("\t", $_, -1) ];
  @$ps >= 18 || die "not 19 fileds:\n$_\n";
  my ($id, $par, $chr, $beg, $end, $srd, $locES, $locIS, $locCS, $loc5S, $loc3S, $phaseS, $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
  $id =~ s/\.TAIR\d+$//;
  $id =~ s/\.v[\d\.]+$//;
  $par =~ s/\.TAIR\d+$//;
  $par =~ s/\.v[\d\.]+$//;
  $ps->[0] = $id;
  $ps->[1] = $par;
  print $fho join("\t", @$ps)."\n";
}
close $fhi;
close $fho;


__END__
