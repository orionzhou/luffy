#!/usr/bin/env perl
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  snp2bed.pl - convert Snp file to Bed file

=head1 SYNOPSIS
  
  snp2bed.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (SNP) file
    -o (--out)    output (BED) file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
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
  open ($fhi, "<$fi") || die "cannot read $fi\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

while( <$fhi> ) {
  chomp;
  next if /(^\#)|(^id)|(^\s*$)/s;
  my ($chr, $pos, $tBase, $qBase, $qid, $qpos, $cid, $lev) = split "\t";
  my $name = "$tBase/$qBase";
  my $score = exists $h_score->{$lev} ? $h_score->{$lev} : 300;
  print $fho join("\t", $chr, $pos-1, $pos, $qBase, $lev)."\n";
}
close $fhi;
close $fho;


__END__
