#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  vcf.anno.pl - run snpEff to annotate a VCF file

=head1 SYNOPSIS
  
  vcf.anno.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input VCF file
    -o (--out)    output Tbl file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Common;
use Gatk;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my $fr = "\$genome/HM101/11_genome.fasta";
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

### java -jar $src/snpEff/snpEff.jar -v mt40 t1.vcf > t2.vcf

my ($fhi, $fho);

if ($fi eq "stdin" || $fi eq "-") {
  $fho = \*STDOUT;
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
  next if(/(^\#)|(^\s*$)/);
  my $line = $_;
  my @ps = split("\t", $line);
  my ($chr, $pos, $ref, $alts, $effs) = @ps[0,1,3,4,7];
  my ($eff, $imp) = ("") x 2;
  if($effs =~ /EFF=(\w+)\((\w+)[\)\|]/) {
    ($eff, $imp) = ($1, $2);
  } else {
    die "$chr:$pos no eff\n";
  }
  print $fho join("\t", $chr, $pos, $eff, $imp)."\n";
}
close $fhi;
close $fho;

__END__
