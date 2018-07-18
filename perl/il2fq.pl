#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  il2fq.pl - convert Fastq-Illumina file to Fastq-Sanger file

=head1 SYNOPSIS
  
  il2fq.pl [-help] [-in input] [-out output]

  Options:
      -h (--help)   brief help message
      -i (--in)     input file
      -o (--out)    output file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my ($fhi, $fho);
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $in = Bio::SeqIO->new(-fh => $fhi, -format => 'fastq-illumina');
my $ou = Bio::SeqIO->new(-fh => $fho, -format => 'fastq');

while (my $seq = $in->next_seq) {
  $ou->write_seq($seq);
}
$in->close();
$ou->close();



exit 0;
