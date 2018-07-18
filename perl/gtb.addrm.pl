#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.addrm.pl - Annotate genes with RepeatMasker output

=head1 SYNOPSIS
  
  gtb.addrm.pl [-help] [-i input-file] [-b overlap-file] [-o output-file]

  Options:
    -h (--help)  brief help message
    -i (--fi)    input file (*.gtb)
    -o (--fo)    output file (*.gtb)
    -b (--fb)    Bed overlap file with RepeatMasker output (*.bed)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gtb;

my ($fi, $fo, $fb) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "fi|i=s"  => \$fi,
  "fo|o=s"  => \$fo,
  "fb|b=s"  => \$fb,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fb;
  
open(my $fhb, $fb) || die "cannot read file $fb\n";
my $h = {};
while(<$fhb>) {
  chomp;
  my ($c1, $b1, $e1, $id, $c2, $b2, $e2) = split "\t";
  $h->{$id} ||= 0;
  $h->{$id} += $e2 - $b2;
}
close $fhb;

open(my $fhi, $fi) || die "cannot read file $fi\n";
open(my $fho, ">$fo") || die "cannot write $fo\n";
print $fho join("\t", @HEAD_GTB)."\n";
my $cnt = 0;
my $h2 = {};
while(<$fhi>) {
  chomp;
  /(^id)|(^\#)|(^\s*$)/i && next;
  my $ps = [ split("\t", $_, -1) ];
  @$ps >= 18 || die "not 19 fileds:\n$_\n";
  my ($id, $par, $chr, $beg, $end, $srd, 
    $locES, $locIS, $locCS, $loc5S, $loc3S, $phase, 
    $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
  my $leng = $end - $beg + 1;
  my $lenr = exists $h->{$id} ? $h->{$id} : 0;
  if($lenr / $leng >= 0.6 && $cat2 ne 'TE') {
    $ps->[15] = 'TE';
    $ps->[16] = 'RepeatMasker';
    $cnt ++;
    $h2->{$cat2} ||= 0;
    $h2->{$cat2} ++;
  }
  print $fho join("\t", @$ps)."\n";
}
close $fhi;
close $fho;
close $fhb;
print "$cnt new TEs added\n";
for my $fam (keys(%$h2)) {
  print "$fam\t$h2->{$fam}\n";
}

