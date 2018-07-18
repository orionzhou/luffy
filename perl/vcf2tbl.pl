#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  vcf2tbl.pl - convert VCF file to TBL file

=head1 SYNOPSIS
  
  vcf2tbl.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (VCF) file
    -o (--out)    output (TBL) file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

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
  next if /(^\#)|(^\s*$)/;
  my ($chr, $pos, $id, $ref, $alt, $qual, $fil, $info, $fmt, @sams) = 
    split "\t";
  my @alts = split(",", $alt);
  my $reflen = length($ref);
  my @altlens = uniq(map {length($_)} @alts);
  @altlens == 1 || die "$chr:$pos $ref-$alt alts not same lens\n";
  my $altlen = $altlens[0];
  
#  ($reflen == 1 || $altlen == 1) || die "$chr:$pos $ref-$alt MNP\n";
  $alt = $alts[0];
  print $fho join("\t", $chr, $pos, $ref, $alt, $qual)."\n";
}
close $fhi;
close $fho;


__END__
