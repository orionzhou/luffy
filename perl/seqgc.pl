#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqgc.pl - calculate GC percentage for a multi-fasta file

=head1 SYNOPSIS
  
  seqgc.pl [-help] [-in input-fasta] [-out output-tbl]

  Options:
    -h (--help)   brief help message
    -i (--in)     input fasta file (can be 'stdin')
    -o (--out)    output (can be 'stdout')

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Time::HiRes qw/gettimeofday tv_interval/;

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
if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');
my $t0 = [gettimeofday];
my $cnt = 0;

while(my $seqO = $seqHI->next_seq()) {
  my ($id, $len, $seq) = ($seqO->id, $seqO->length, $seqO->seq);
  print $fho join("\t", $id, calc_gc($seq))."\n";
  
  $cnt ++;
  printf "%d: %.01f min\n", $cnt, tv_interval($t0, [gettimeofday]) / 60 if $cnt % 100000 == 0;
}
close $fhi;
close $fho;

sub calc_gc {
  my ($str) = @_;
  my ($cntGC, $cntN) = (0, 0);
  while($str =~ /([GCN])/ig) {
    if(uc($1) eq "N") {
      $cntN ++;
    } else {
      $cntGC ++;
    }
  }
  my $len = length($str) - $cntN;
  return 0 if $len == 0;
  return sprintf "%.03f", $cntGC / $len;
}

exit 0;
