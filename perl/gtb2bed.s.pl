#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2bed.pl - convert a Gtb file to BED format

=head1 SYNOPSIS
  
  gtb2bed.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gtb format)
    -o (--out)    output file (Bed format)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Common;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

  open(my $fhi, $fi) || die "cannot read file $fi\n";
  open(my $fho, ">$fo") || die "cannot write $fo\n";
  while(<$fhi>) {
    chomp;
    /(^id)|(^\#)|(^\s*$)/i && next;
    my $ps = [ split("\t", $_, -1) ];
    @$ps >= 18 || die "not 19 fileds:\n$_\n";
    my ($id, $par, $chr, $beg, $end, $srd, 
      $locES, $locIS, $locCS, $loc5S, $loc3S, $phase, 
      $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
    print $fho join("\t", $chr, $beg-1, $end, $id)."\n";
  }
close $fhi;
close $fho;

