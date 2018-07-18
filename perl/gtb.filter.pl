#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.filter.pl - filter a Gtb file 

=head1 SYNOPSIS
  
  gtb.filter.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (Gtb) file
    -o (--out)    output (Gtb) file
    -l (--len)    minimum CDS length (default: 30)
    -c (--chr)    file with Chr IDs

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

my ($fi, $fo, $fc) = ('') x 3;
my $min_len = 30;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "len|l=i"  => \$min_len,
  "chr|c=s"  => \$fc,
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

my $hc;
if($fc) {
  open(my $fhc, "<$fc") or die "cannot read $fc\n";
  while(<$fhc>) {
    chomp;
    /(^\#)|(^\s*$)/ && next;
    $hc->{$_} = 1;
  }
  close $fhc;
}

my ($cntf, $cnta) = (0, 0);
print $fho join("\t", @HEAD_GTB)."\n";
while( <$fhi> ) {
  chomp;
  /(^id)|(^\#)|(^\s*$)/i && next;
  my $ps = [ split("\t", $_, -1) ];
  @$ps >= 18 || die "not 19 fileds:\n$_\n";
  my ($id, $par, $chr, $beg, $end, $srd, $locES, $locIS, $locCS, $loc5S, $loc3S, $phaseS, $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
  my $flag = 1;

  $flag = 0 if $conf < $min_len;
  $flag = 0 if $fc && exists $hc->{$chr};

  $cnta ++;
  $cntf ++ if $flag == 0;
  print $fho join("\t", @$ps)."\n" if $flag == 1;
}
close $fhi;
close $fho;

printf "%5d | %5d filtered\n", $cntf, $cnta;

__END__
