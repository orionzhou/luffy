#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.addpfam.pl - add PFAM annotation to GTB 

=head1 SYNOPSIS
  
  gtb.addpfam.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (*.gtb)
    -o (--out)    output file (*.gtb)
    -p (--pfam)   PFAM hmmscan output (*.tbl)
    -c (--cat)    genefam category file (*.tbl)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Gtb;
use List::Util qw/min max sum/;

my ($fi, $fo, $fp) = ('') x 3;
my $fc = "$ENV{'data'}/db/pfam/genefam.tsv";
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "pfam|p=s" => \$fp,
  "cat|c=s"  => \$fc,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fp;

my ($fhi, $fho);
if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot read $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

open(my $fhc, "<$fc") or die "cannot read $fc\n";
my $hd;
while(<$fhc>) {
  chomp;
  /^fam/i && next;
  my ($fam, $dom) = split "\t";
  $hd->{$dom} = $fam;
}
close $fhc;

open(my $fhp, "<$fp") or die "cannot read $fp\n";
my $h;
while(<$fhp>) {
  chomp;
  /^qid/i && next;
  my @ps = split "\t";
  my ($qid, $hid, $e) = @ps[0, 4, 8];
  $h->{$qid}->{$hid} ||= $e;
  $h->{$qid}->{$hid} = min($e, $h->{$qid}->{$hid});
}
close $fhp;

my $hs;
my $cnt = 0;
print $fho join("\t", @HEAD_GTB)."\n";
while( <$fhi> ) {
  chomp;
  /(^id)|(^\#)|(^\s*$)/i && next;
  my $ps = [ split("\t", $_, -1) ];
  @$ps >= 18 || die "not 19 fileds:\n$_\n";
  my ($id, $par, $chr, $beg, $end, $srd, $locES, $locIS, $locCS, $loc5S, $loc3S, $phaseS, $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
  if(exists $h->{$id}) {
    $cnt ++;
    my $h2 = $h->{$id};
    my @doms = sort {$h2->{$a} <=> $h2->{$b}} keys(%$h2);
    $cat3 = $doms[0];
    $cat2 = exists $hd->{$doms[0]} ? $hd->{$doms[0]} : 'Miscellaneous';
    $note = join(" ", @doms);
  } else {
    $cat2 = 'Unknown';
    $note = '';
  }
  $ps->[15] = $cat2;
  $ps->[16] = $cat3;
  $ps->[17] = $note;
  $hs->{$cat2} ||= 0;
  $hs->{$cat2} ++;
  print $fho join("\t", @$ps)."\n";
}
close $fhi;
close $fho;
printf "%6d models annotated with Pfam info\n", $cnt;
for my $cat (sort {$hs->{$b} <=> $hs->{$a}} keys(%$hs)) {
  printf "%6d: %s\n", $hs->{$cat}, $cat;
}


__END__
