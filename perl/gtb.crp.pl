#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.crp.pl - 

=head1 SYNOPSIS
  
  gtb2bed.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gtb format)
    -o (--out)    output file (Gtb format)

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
  my $ti = readTable(-in => $fi, -header => 1);
  $ti->delCol("seq");
  for my $i (0..$ti->lastRow) {
    my ($id, $par, $chr, $beg, $end, $srd, $locE, $locI, $locC, $loc5, $loc3, $phase, $src, $conf, $cat1, $cat2, $cat3, $note) = $ti->row($i);
    $ti->setElm($i, "cat1", "mRNA");
    $ti->setElm($i, "note", "");
    if($cat3 ge "CRP0000" && $cat3 le "CRP1030") {
      $cat2 = "CRP0000-1030";
    } elsif($cat3 ge "CRP1040" && $cat3 le "CRP1530") {
      $cat2 = "NCR";
    } elsif($cat3 ge "CRP1600" && $cat3 le "CRP6250") {
      $cat2 = "CRP1600-6250";
    } else {
      die "unknonw cat: $cat3\n";
    }
    $ti->setElm($i, "cat2", $cat2);
  }
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $ti->tsv(1);
  close $fho;


