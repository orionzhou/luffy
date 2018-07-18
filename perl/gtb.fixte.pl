#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.fixte.pl - fix TE annotation

=head1 SYNOPSIS
  
  gtb.fixte.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (Gtb) file
    -o (--out)    output (Gtb) file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Gtb;
use Pod::Usage;
use Common;

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
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open($fhi, "<$fi") or die "cannot read $fi\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open($fho, ">$fo") or die "cannot write $fo\n";
}

print $fho join("\t", @HEAD_GTB)."\n";
while( <$fhi> ) {
  chomp;
  /(^id)|(^\#)|(^\s*$)/i && next;
  my $ps = [ split("\t", $_, -1) ];
  @$ps >= 18 || die "not 19 fileds:\n$_\n";
  my ($id, $par, $chr, $beg, $end, $srd, $locES, $locIS, $locCS, $loc5S, $loc3S, $phaseS, $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
  if($cat1 =~ /^gene$/i) {
    $ps->[16] = "gene";
  } elsif($cat1 =~ /transposable_element/i) {
    $ps->[16] = "TE";
  }
  print $fho join("\t", @$ps)."\n";
}
close $fhi;
close $fho;



__END__
