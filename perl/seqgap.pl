#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqgap.pl - report gap ('N's) locations in fasta sequences

=head1 SYNOPSIS
  
  seqgap.pl [-help] [-min minimum-gap] [-in input-file] [-out output-file] 

  Options:
    -h (--help)   brief help message
    -i (--in)     input (FASTA format) (default: 'stdin')
    -o (--out)    output (BED format) (default: "stdout') 
    -m (--min)    minimum length of gaps to report (default: 10)

=head1 DESCRIPTION
  
  For each gap (consecutive 'N's) in input sequences, report its locations
  in tabular-separated  BED format with three fields: 'sequence-ID', 
  'start' and 'end'.
  Example output:
    chr1   5      50
    chr1   100    131
    chr5   56     100

=head1 AUTHOR
  
  Peng Zhou - University of Minnesota

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($fi, $fo) = ('') x 2;
my $fhi;
my $fho;
my $len_min = 10;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "min|m=i" => \$len_min,
) or pod2usage(2);
pod2usage(1) if $help_flag;

# read from stdin or opens an input file handler
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

# write to stdout or opens an output file handler
if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

# initialize a Bio::SeqIO object using input handler
my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');

# iterate through input sequences
while(my $seqO = $seqHI->next_seq()) {
  my ($id, $seq, $seqlen) = ($seqO->id, $seqO->seq, $seqO->length);

# using regexp to find gaps (consecutive 'N's in a string),
# determine their locations and lengths, filter, and output
  while($seq =~ /N+/ig) {
    my ($beg, $end) = ($-[0]+1, $+[0]);
    my $len = $end - $beg + 1;

# don't report short gaps
    $len >= $len_min || next;

# note that BED format use 0-based start coordinates
    print $fho join("\t", $id, $beg - 1, $end)."\n";
  }
}
$seqHI->close();
close $fho;

exit 0;
