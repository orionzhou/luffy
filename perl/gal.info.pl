#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal.info.pl - brief summary

=head1 SYNOPSIS
  
  gal.info.pl [-help] [-in input-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file

=cut
  
#### END of POD documentation.
#--------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fi, $fo) = ('') x 2;
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi;

if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

my $h;
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps == 21;
  my ($qId) = @$ps[1];
  $h->{$qId} ||= 0;
  $h->{$qId} ++;
}
close $fhi;

my $cntQ = scalar(keys %$h);
my $cntT = sum(values(%$h));
printf "%8d queries -> %8d targets\n", $cntQ, $cntT;
my @dups = grep {$h->{$_} > 1} keys %$h;
my $cntDupHits = @dups ? sum( map {$h->{$_}} @dups ) : 0;
my $cntU = $cntQ - @dups;
printf "    %8d uniquely mapped\n", $cntU;
printf "    %8d non-uniquely mapped to %8d targets\n", scalar(@dups), $cntDupHits;


__END__
