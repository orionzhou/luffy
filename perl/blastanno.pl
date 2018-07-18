#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  blastanno.pl - annotate BLAST output using NCBI taxonomy database

=head1 SYNOPSIS
  
  blastanno.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file

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
use Eutils;

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

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $t = readTable(-in => $fi, -header => 1);
my @gis;
for my $i (0..$t->nofRow-1) {
  my $id = $t->elm($i, "tId");
  if($id =~ /gi\|(\d+)\|/) {
    push @gis, $1;
    $t->setElm($i, "tId", $1);
  } else {
    die "unknown tId: $id\n";
  }
}
my $h1 = gi2Taxid(@gis);
my $h2 = annotate_taxid(values(%$h1));

$t->addCol([('')x$t->nofRow], 'superkingdom');
$t->addCol([('')x$t->nofRow], 'kingdom');
$t->addCol([('')x$t->nofRow], 'family');
$t->addCol([('')x$t->nofRow], 'species');
$t->addCol([('')x$t->nofRow], 'cat');
for my $i (0..$t->nofRow-1) {
  my $gi = $t->elm($i, "tId");
  my $taxid = $h1->{$gi};
  my @cats = @{$h2->{$taxid}};
  my $cat = $cats[1] eq "Viridiplantae" ? "plant" : "foreign";
  $t->setElm($i, "superkingdom", $cats[0]);
  $t->setElm($i, "kingdom", $cats[1]);
  $t->setElm($i, "family", $cats[2]);
  $t->setElm($i, "species", $cats[3]);
  $t->setElm($i, "cat", $cat);
}

print $fho $t->tsv(1);
close $fho;

