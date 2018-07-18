#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  net2gal.pl - convert a Net/Chain file to Gal format

=head1 SYNOPSIS
  
  net2gal.pl [-help] [-net net-file] [-chain chain-file] [-out output]

  Options:
    -h (--help)   brief help message
    -n (--net)    Net file
    -c (--chain)  chain file (netChainSubset result)
    -o (--out)    output (Gal) file

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

my ($fn, $fc, $fo) = ('') x 3;
my ($fhn, $fhc, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"     => \$help_flag,
  "net|n=s"    => \$fn,
  "chain|c=s"  => \$fc,
  "out|o=s"    => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fn || !$fc || !$fo;

if ($fc eq "stdin" || $fc eq "-") {
  $fhc = \*STDIN;
} else {
  open ($fhc, $fc) || die "Can't read file $fc: $!\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $h = {};
open ($fhn, $fn) || die "Can't read file $fn\n";
my $hi;
my $ph = { 0 => "" };
while( <$fhn> ) {
  chomp;
  next if /(^\#)|(^\s*$)/;
  if(/(^ +)(fill.*)$/) {
    my $lev = (length($1)+1)/2;
    my @ps = split(" ", $2);
    my $n = (@ps-7) / 2;
    my %h = map {$ps[7+$_*2] => $ps[7+$_*2+1]} (0..$n-1);
    my $id = $h{"id"};
    my $type = $h{"type"};
    my $qDup = exists $h{"qDup"} ? $h{'qDup'} : '';
    my $qFar = exists $h{'qFar'} ? $h{'qFar'} : '';
    my $qOver = exists $h{'qOver'} ? $h{'qOver'} : '';
    
    exists $ph->{$lev-1} || die "no level $lev-1\n";
    my $pid = $ph->{$lev-1};

    $hi->{$id} ||= 0;
    $hi->{$id} ++;
    $id .= ".".$hi->{$id} if $hi->{$id} > 1;

    $h->{$id} = [$lev, $type, $pid, $qDup, $qFar, $qOver];

    $ph->{$lev} = $id;
  }
}
close $fhn;


$hi = {};
print $fho join("\t", @HEAD_GAL)."\n";
while( <$fhc> ) {
  chomp;
  my @ps = split /\s/;
  if($ps[0] eq "chain") {
    my ($score, $tId, $tSize, $tSrd, $tBeg, $tEnd, 
      $qId, $qSize, $qSrd, $qBeg, $qEnd, $id) = @ps[1..$#ps];
    $tSrd eq "+" || die "$id: tSrd -\n";
    $tBeg += 1;
    $qBeg += 1;
    ($qBeg, $qEnd) = ($qSize-$qEnd+1, $qSize-$qBeg+1) if $qSrd eq "-";
    
    $hi->{$id} ||= 0;
    $hi->{$id} ++;
    $id .= ".".$hi->{$id} if $hi->{$id} > 1;
    
    my ($td, $qd) = (0, 0);
    my ($rtloc, $rqloc) = ([], []);
    while( <$fhc> ) {
      last if /^\s*\n$/;
      my @pps = split /\s/;
      my $len = $pps[0];
      my ($dt, $dq) = (0, 0);
      ($dt, $dq) = @pps[1..2] if @pps >= 3;
      push @$rtloc, [$td+1, $td+$len];
      push @$rqloc, [$qd+1, $qd+$len];
      
      $td += $len + $dt;
      $qd += $len + $dq;
    }
    my ($rtlocs, $rqlocs) = (locAry2Str($rtloc), locAry2Str($rqloc));

    exists $h->{$id} || die "no stat for $id\n";
    my ($lev, $type, $pid, $qDup, $qFar, $qOver) = @{$h->{$id}};
    print $fho join("\t", $id, $qId, $qBeg, $qEnd, $qSrd, $qSize,
      $tId, $tBeg, $tEnd, $tSrd, $tSize, 
      ('') x 6, $score, 
      $lev, $type, $pid, $qDup, $qFar, $qOver, $rqlocs, $rtlocs)."\n";
  }
}
close $fhc;
close $fho;


__END__
