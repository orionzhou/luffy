#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal.addlev.pl - fill in level field

=head1 SYNOPSIS
  
  gal.addlev.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gal)
    -o (--out)    output file (Gal)
    -n (--net)    net file

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

my ($fi, $fo, $fn) = ('') x 3;
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "net|n=s"  => \$fn,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fn;

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

open(my $fhn, "<$fn") || die "cannot read $fn\n";
my ($hl, $hi);
my $ph = {0 => ""};
while( <$fhn> ) {
  chomp;
  next if /(^\#)|(^\s*$)/;
  if(/(^ +)(fill.*)$/) {
    my $lev = (length($1) + 1) / 2;
    my @ps = split(" ", $2);
    my $n = (@ps - 7) / 2;
    my %h = map {$ps[7+$_*2] => $ps[7+$_*2+1]} (0..$n-1);
    my $id = $h{"id"};
    my $type = $h{"type"};
    exists $ph->{$lev-1} || die "no level $lev-1\n";
    my $pid = $ph->{$lev-1};

    $hi->{$id} ||= 0;
    $hi->{$id} ++;
    $id .= ".".$hi->{$id} if $hi->{$id} > 1;

    $hl->{$id} = $lev;
    $ph->{$lev} = $id;
  }
}
close $fhn;

print $fho join("\t", @HEAD_GAL)."\n";
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my @ps = split "\t";
  next unless @ps == 21;
  my ($cid, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $lev, $ali, $mat, $mis, $qN, $tN, $ident, $score, $tlS, $qlS) = @ps;
  exists $hl->{$cid} || die "no lev for $cid\n";
  $ps[11] = $hl->{$cid};
  print $fho join("\t", @ps)."\n";
}
close $fhi;
close $fho;


__END__
