#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 SYNOPSIS
  
  winsplit.pl [-help] [-in window-file] [-n number-chunks] [-out output]

  Options:
    -help   brief help message
    -in     window file (can be 'stdin')
    -out    output (can be 'stdout') 
    -n      number of chunks (default: 1)
    -part   part number (default: 1)
    -ovlp   overlap btw. windows (default: 0)

=head1 BUGS
  
=head1 REFERENCES
  
=head1 VERSION
  
  0.1
  
=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use POSIX qw/ceil floor/;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my ($n, $part, $ovlp) = (1, 1, 0);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "n=i"      => \$n,
  "part|p=i" => \$part,
  "ovlp|v=i" => \$ovlp,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if $n < 1 || $part < 1 || $part > $n || $ovlp < 0;

my ($fhi, $fho);
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

my @stats;
my $lentot = 0;
while(<$fhi>) {
  chomp;
  next if /(^#)|(^id\s)|(^chr\s)/;
  my @ps = split "\t";
  my ($chr, $beg, $end) = ($ps[0], '', '');
  next unless @ps >= 2;
  if(@ps == 2) {
    ($beg, $end) = (1, $ps[1]);
  } else {
    ($beg, $end) = @ps[1,2];
  }
  my $len = $end - $beg + 1;
  push @stats, [$chr, $beg, $end, $len, $lentot+$beg, $lentot+$end];
  $lentot += $len;
}
close $fhi;

my $size = ceil( ($lentot - $ovlp*($n+1)) / $n );

#print $lentot."\t".$size."\n";

my $bega = ($ovlp+$size) * ($part-1) + 1;
my $enda = ($ovlp+$size) * $part + $ovlp;
$enda = $lentot if $part == $n;
#print $fho join("\t", $part, $bega, $enda, $enda-$bega+1)."\n";

for (@stats) {
  my ($chr, $begir, $endir, $len, $begia, $endia) = @$_;
  if( max($begia,$bega) <= min($endia, $enda) ) {
    my $begoa = max($begia, $bega);
    my $endoa = min($endia, $enda);
    my $begor = $begir + ($begoa - $begia);
    my $endor = $begir + ($endoa - $begia);
    print $fho join("\t", $chr, $begor, $endor)."\n";
  }
}

exit 0;



