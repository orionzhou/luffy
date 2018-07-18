#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  idm2bed.pl - convert Idm file to EqualIndel file

=head1 SYNOPSIS
  
  idm2bed.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input (IDM) file
      -out    output file
      -pre    ID prefix

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;

my ($fi, $fo) = ('') x 2;
my $pre = "pre";
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "pre|p=s"  => \$pre,
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
print $fho join("\t", qw/name chromosome start end allele ref_allele/)."\n";

my $i = 0;
while( <$fhi> ) {
  chomp;
  next if /(^\#)|(^\s*$)/;
  my ($chr, $beg, $end, $tBase, $qBase) = split "\t";
  my $tlen = length($tBase);
  my $qlen = length($qBase);
  $end - $beg - 1 == $tlen || die "len error: $chr:$beg-$end\[$tBase]\n";
  my $id = $pre . ++$i;
  $tBase = "-" if $tlen == 0;
  $qBase = "-" if $qlen == 0;

  print $fho join("\t", $id, $chr, $beg + 1, $end - 1, $qBase, $tBase)."\n";
}
close $fhi;
close $fho;


__END__
