#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqconcat.pl - concatenate all sequences from input file to one sequence

=head1 SYNOPSIS
  
  seqconcat.pl [-help] [-id seq ID] [-gap gap-length] [-in input-file] [-out output-file] [-path path-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input sequence file (fasta)
    -o (--out)    output sequence file (fasta)
    -p (--path)   output path file (tabular)
    -d (--id)     output sequene ID (default: "chrU")
    -g (--gap)    gap length (default: 1000)

=cut
  
#### END of POD documentation.
#------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($fi, $fo, $fp) = ('') x 3;
my $gapLen = 1000;
my $seqid = "chrU";
my $help_flag;

#------------------------------ MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "path|p=s" => \$fp,
  "gap|g=i" => \$gapLen,
  "id|d=s"  => \$seqid,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fp;

open(my $fhi, "<$fi") || die "cannot read $fi";
open(my $fho, ">$fo") || die "cannot write $fo";
open(my $fhp, ">$fp") || die "cannot write $fp";
print $fhp join("\t", qw/id chr beg end strand/)."\n";

my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');
my ($endPrev, $cnt) = (0, 0);
my $seqStr = "";
while(my $seqO = $seqHI->next_seq()) {
  my ($seqLen, $seq) = ($seqO->length, $seqO->seq);
  if($endPrev > 0) {
    $seqStr .= 'N' x $gapLen;
    $endPrev += $gapLen;
  }
  print $fhp join("\t", $seqO->id, $seqid, $endPrev+1, $endPrev+$seqLen, "+")."\n";
  $seqStr .= $seq;
  $endPrev += $seqLen;
}
$seqHI->close();

my $seqHO = Bio::SeqIO->new(-fh=>$fho, -format=>'fasta');
$seqHO->write_seq(Bio::Seq->new(-id=>$seqid, -seq=>$seqStr));
$seqHO->close();

exit 0;


