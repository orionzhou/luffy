#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  chain2gal.pl - convert a Chain file to Gal format

=head1 SYNOPSIS
  
  chain2gal.pl [-help] [-in input-file] [-out output-file]

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

my $hi;
print $fho join("\t", @HEAD_GAL)."\n";
while( <$fhi> ) {
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
    my $ali = 0;
    my ($rtloc, $rqloc) = ([], []);
    while( <$fhi> ) {
      last if /^\s*\n$/;
      my @pps = split /\s/;
      my $len = $pps[0];
      my ($dt, $dq) = (0, 0);
      ($dt, $dq) = @pps[1..2] if @pps >= 3;
      push @$rtloc, [$td+1, $td+$len];
      push @$rqloc, [$qd+1, $qd+$len];
     
      $ali += $len;
      $td += $len + $dt;
      $qd += $len + $dq;
    }
    my ($rtlocs, $rqlocs) = (locAry2Str($rtloc), locAry2Str($rqloc));
    print $fho join("\t", $id, $tId, $tBeg, $tEnd, $tSrd, $tSize,
      $qId, $qBeg, $qEnd, $qSrd, $qSize, 
      '', $ali, $ali, 0, 0, 0, '', $score, $rtlocs, $rqlocs)."\n";;
  }
}
close $fhi;
close $fho;


__END__
