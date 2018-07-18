#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  sammapp.pl - generate report of mapping uniqueness for a SAM file

=head1 SYNOPSIS
  
  sammapp.pl [-help] [-in input] [-out output]

  Options:
      -help   brief help message
      -in     input SAM file (can be 'stdin')
      -out    output file (can be 'stdout')
      -mis    maximum tolerance of mapping distance (default: 3)

=head1 BUGS
  
=head1 REFERENCES
  
=head1 VERSION
  
  0.1
  
=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my $nm_max = 3;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "mismatch|m=i" => \$nm_max,
) or pod2usage(2);
pod2usage(1) if $help_flag;
my ($fhi, $fho);

if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $pid = '';
my $h = { map {$_ => 0} (0..$nm_max) };
while( <$fhi> ) {
  chomp;
  next if /(^\#)|(^\@)/;
  my @ps = split "\t";
  my $id = $ps[0];
  my ($nm, $yf) = ('') x 2;
  for my $i (11..$#ps) {
    if($ps[$i] =~ /NM\:i\:(\d+)/) {
      $nm = $1;
      last;
    }
}
#    die "no edit distance for $id [YF:Z:$yf]\n" if ($nm eq '' && $yf ne "NS");

  if($pid eq '' || $pid eq $id) {
    $pid = $id if $pid eq '';
    $h->{$nm} ++ if exists $h->{$nm};
  } else {
    my $n = sum( map {$h->{$_}} (0..$nm_max) );
    my $mapp = $n == 0 ? 0 : sprintf("%.03f", 1/$n);
    print $fho join("\t", $pid, $mapp)."\n";
    $h = { map {$_ => 0} (0..$nm_max) };

    $pid = $id;
    $h->{$nm} ++ if exists $h->{$nm};
  }
}
my $n = sum( map {$h->{$_}} (0..$nm_max) );
my $mapp = $n == 0 ? 0 : sprintf("%.03f", 1/$n);
print $fho join("\t", $pid, $mapp)."\n";

close $fhi;
close $fho;

exit 0;



