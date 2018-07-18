#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  vcf.fix.indel.pl - LeftAlignAndTrim InDels

=head1 SYNOPSIS
  
  vcf.fix.indel.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input VCF file
    -o (--out)    output file
    -r (--ref)    ref-fasta file (default: $genome/HM101/11_genome.fasta)

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

my ($fhi, $fho);
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
  if(/(^\#)|(^\s*$)/) {
    print $fho $_."\n";
    next;
  }
  my $line = $_;
  my @ps = split("\t", $line);
  my ($ref, $alts) = @ps[3,4];
  my @alts = split(",", $alts);
  if(@alts == 1) {
    print $fho $line."\n";
    next;
  }
  @alts == 2 || die ">2 alts: $line\n";
  my ($alt1, $alt2) = @alts;
  ($alt1, $alt2) = ($alt2, $alt1) if length($alt1) > length($alt2);
  my ($lenr, $lena1, $lena2) = map {length($_)} ($ref, $alt1, $alt2);
  my ($cr, $ca1, $ca2) = map {substr($_, 0, 1)} ($ref, $alt1, $alt2);
  my $alt;
  if($lenr == 1 && $lena1 == 1 && $lena2 > 1) {
    die "not ins+snp: $ref $alt1 $alt2\n$line\n" if $ref eq $alt1 || $ref ne substr($alt2, 0, 1);
    $alt = $alt1 . substr($alt2, 1);
  } elsif($lenr > 1 && $lenr == $lena2 && $lena1 == 1) {
    die "not del+snp: $ref $alt1 $alt2\n$line\n" if substr($ref, 0, 1) ne $alt1 || substr($ref, 0, 1) eq substr($alt2, 0, 1) || substr($ref, 1) ne substr($alt2, 1);
    $alt = substr($alt2, 0, 1);
  } elsif($lenr > 1 && $lenr < $lena2 && $lena1 == 1) {
    print join("\t", $ref, $alt1, $alt2)."\n";
    $alt = $alt2;
  } elsif($lenr == 1 && (($ref eq $ca1) || ($ref eq $ca2))) {
    print join("\t", $ref, $alt1, $alt2)."\n";
  } else {
    die "error: $line\n";
  }
  $ps[4] = $alt;
  $ps[10] = "1/1";
  print $fho join("\t", @ps)."\n";
}
close $fhi;
close $fho;

#runCmd("java -jar $gatk -T LeftAlignAndTrimVariants \\
#  -R $fr --trimAlleles \\
#  --variant $fi -o $fo");
#runCmd("rm $fi*.idx $fo*.idx");


__END__
