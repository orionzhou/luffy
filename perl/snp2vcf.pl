#!/usr/bin/env perl
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  snp2vcf.pl - convert Snp file to Vcf file

=head1 SYNOPSIS
  
  snp2vcf.pl [-help] [-in input-file] [-out output-file]

  Options:
    -help   brief help message
    -in     input (SNP) file
    -out    output (VCF) file
    -sample sample ID 

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Vcfhead;

my ($fi, $fo) = ('') x 2;
my $sample = "";
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "sample|s=s"  => \$sample,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$sample;

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

print $fho $vcfhead."\n";
print $fho "#".join("\t", @colhead, $sample)."\n";

while( <$fhi> ) {
  chomp;
  next if /(^\#)|(^id)|(^\s*$)/s;
  my ($chr, $pos, $tBase, $qBase, $qid, $qpos, $cid, $lev) = split "\t";
  $lev == 1 || next;
  print $fho join("\t", $chr, $pos, ".", $tBase, $qBase, 50, '.',
    '.', 'GT', '1/1')."\n";
}
close $fhi;
close $fho;
#runCmd("bgzip -f $fo");
#runCmd("tabix -p vcf $fo.gz");


__END__
