#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtbx2fas.pl - convert a Gtbx file to fasta format

=head1 SYNOPSIS
  
  gtbx2fas.pl [-help] [-in input-file] [-out output-file]

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
use Common;
use Location;
use Bio::Seq;
use Bio::SeqIO;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my ($fhi, $fho);
if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot open $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $seqH = Bio::SeqIO->new(-fh => $fho, -format => "fasta");
my $t = readTable(-inh => $fhi, -header => 1);

for my $i (0..$t->lastRow) {
  my ($id, $seqstr) = map {$t->elm($i, $_)} qw/id seq/;
  my $seq = Bio::Seq->new(-id => $id, -seq => $seqstr); 
  $seqH->write_seq($seq);
}
$seqH->close();

__END__
