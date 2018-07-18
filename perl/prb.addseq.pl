#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  prb.addseq.pl - add sequance for a probe file

=head1 SYNOPSIS
  
  prb.addseq.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (TBL)
    -o (--out)    output file (TBL)
    -s (--seq)    sequence db (default: $genome/pan4/11_genoeme.fas)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::Fasta;
use List::Util qw/min max sum/;

#--------------------------------- MAIN -----------------------------------#
my ($fi, $fo) = ('', '');
my $fs = "/home/youngn/zhoup/Data/genome/pan4/11_genome.fas";
my $help_flag;
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "seq|s=s" => \$fs,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fs;

my ($fhi, $fho);
if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $db = Bio::DB::Fasta->new($fs); 

while(<$fhi>) {
  chomp;
  my @ps = split "\t";
  if($ps[0] =~ /^chr/i) {
    print $fho join("\t", @ps, 'seq')."\n";
    next;
  }
  my ($chr, $beg, $end) = @ps;
  my $seq = $db->seq($chr, $beg, $end);
  print $fho join("\t", @ps, $seq)."\n";
}
close $fhi;
close $fho;

exit 0;
