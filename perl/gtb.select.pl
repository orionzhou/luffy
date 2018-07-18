#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.filter.pl - filter a Gtb file 

=head1 SYNOPSIS
  
  gtb.select.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (Gtb) file
    -o (--out)    output (Gtb) file
    -c (--chr)    Chr IDs separated by comma (,)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Gtb;
use Common;
use Location;

my ($fi, $fo, $cids) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "chr|c=s"  => \$cids,
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

my @cids = split(",", $cids);
my $hc = { map {$_ => 1} @cids };
print join(" ", keys(%$hc))."\n";

my ($cnts, $cnta) = (0, 0);
print $fho join("\t", @HEAD_GTB)."\n";
while( <$fhi> ) {
  chomp;
  /(^id)|(^\#)|(^\s*$)/i && next;
  my $ps = [ split("\t", $_, -1) ];
  @$ps >= 18 || die "not 19 fileds:\n$_\n";
  my ($id, $par, $chr, $beg, $end, $srd, $locES, $locIS, $locCS, $loc5S, $loc3S, $phaseS, $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
  my $flag = 0;
  $flag = 1 if exists $hc->{$chr};

  $cnta ++;
  $cnts ++ if $flag == 1;
  print $fho join("\t", @$ps)."\n" if $flag == 1;
}
close $fhi;
close $fho;

printf "%5d | %5d selected\n", $cnts, $cnta;

__END__
