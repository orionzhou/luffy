#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  bedseq.pl - retrieve sequences with given intervals from a seq-db

=head1 SYNOPSIS
  
  bedseq.pl [-help] [-db fasta-db] [-in input-BED] [-out output-BED]

  Options:
    -h (--help)   brief help message
    -d (--db)     sequence database (fasta)
    -i (--in)     input BED file with genomic intervals
    -o (--out)    output file (BED)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Bio::DB::Fasta;
use Bio::SeqIO;

my ($fd, $fo, $fi) = ('') x 3;
my ($fho);
my $help_flag;

#------------------------------ MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "db|d=s"   => \$fd,
  "out|o=s"  => \$fo,
  "in|i=s"   => \$fi,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fd || !$fi || !$fo;

my $db = Bio::DB::Fasta->new($fd);

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $cnt = 0;

open(my $fhi, "<$fi") or die "cannot read $fi\n";
while(<$fhi>) {
  chomp;
  my @ps = split "\t";
  my ($seqid, $beg, $end) = @ps[0..2];
  $beg += 1;
  $beg <= $end || die "loc error in $fi\n$seqid:$beg-$end\n";
  
  my $seq = $db->seq($seqid, $beg, $end);
  defined $seq || die "$seqid not in db\n";
  print $fho join("\t", @ps, $seq)."\n";
  $cnt ++;
}
close $fhi;

printf "  %4d sequences extracted\n", $cnt;

exit 0;



