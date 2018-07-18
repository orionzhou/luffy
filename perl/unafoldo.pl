#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  unafoldo.pl - process unafold(hybrid-ss-min) output

=head1 SYNOPSIS
  
  unafoldo.pl [-help] [-in input-file (*.dG)] [-seq fasta-file] [-out output]

  Options:
    -help   brief help message
    -in     input file (*.dG)
    -seq    fasta file
    -out    output

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($fi, $fs, $fo) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "seq|s=s" => \$fs,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fs;

my ($fhi, $fho);
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi: $!\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo: $!\n";
}

my $seqHI = Bio::SeqIO->new(-file=>"<$fs", -format=>'fasta');
while(<$fhi>) {
  chomp;
  next if /^\#/;
  my @ps = split "\t";
  $ps[0] == 37 || die "error content: ".join("\t", @ps)."\n";
  my $id = $seqHI->next_seq()->id;
  print $fho join("\t", $id, $ps[1])."\n";
}
$seqHI->close();
close $fhi;
close $fho;

exit 0;
