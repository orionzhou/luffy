#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  snp2phy.pl - convert a snp file to phylip file

=head1 SYNOPSIS
  
  snp2phy.pl [-help] [options] [-in input] [-out output] 

  Options:
    -h (--help)  brief help message
    -i (--in)    input file (SNP)
    -o (--out)   output file( Phylip)
    -r (--ref)   whether include reference (ID) in the output 

=head1 VERSION
  
  0.1
  
=cut
  
#### END of POD documentation.
#--------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;

my ($fi, $fo) = ('') x 2;
my $refname = "";
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "ref|r=s" => \$refname,
) or pod2usage(2);
pod2usage(1) if $help_flag;

if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot read $fi\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my (@names, @poss, @data);
my ($nind, $npos) = (0, 0);
my $refseq = "";

my $first = 1;
while( <$fhi> ) {
  chomp;
  my $line = $_;
  my @ps = split("\t", $line);
  if($first) {
    $nind = @ps - 4;
    @data = map {""} (0..$nind-1);
  } else {
    die "not $nind + 4 cols\n" unless $nind + 4 == @ps;
  }
  my ($ref, $alt) = ($ps[2], $ps[3]);
  die "ref[$ref] alt[$alt] not SNP\n" if length($ref) != 1 || length($alt) != 1;
  push @poss, $ps[1];
  $npos ++;

  $refseq .= $ref;
  for my $i (0..$nind-1) {
    if($ps[$i+4] =~ /([A-Za-z0-9\-_]+)\=([01\.])\/([01\.])/) {
      push @names, $1 if $first;
      my $nt = ($2 eq "0" && $3 eq "0") ? $ref : ($2 eq "1" && $3 eq "1") ? $alt : "N";
      $data[$i] .= $nt;
    } else {
      die "unknown allele: $ps[$i+4]\n";
    }
  }
  $first = 0 if $first;
}
close $fhi;

if($refname) {
  print $fho join(" ", $nind + 1, $npos)."\n";
} else {
  print $fho join(" ", $nind, $npos)."\n";
}
for my $r (0..ceil($npos/250)-1) {
  my $refstr = substr($refseq, $r*250, 250);
  if($refname) {
    if($r == 0) {
      printf $fho "%-20s%s\n", $refname, $refstr;
    } else {
      printf $fho "%s\n", $refstr;
    }
  }
  for my $i (0..$nind-1) {
    my $str = substr($data[$i], $r*250, 250);
    if($r == 0) {
      printf $fho "%-20s%s\n", $names[$i], $str;
    } else {
      printf $fho "%s\n", $str;
    }
  }
  print $fho "\n";
}
close $fho;



