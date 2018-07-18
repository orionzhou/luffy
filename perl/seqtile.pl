#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqtile.pl - generate sequcne tile-ups for an input fasta file

=head1 SYNOPSIS
  
  seqtile.pl [-help] [-in input] [-step win-step] [-size win-size] [-out output]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (can be 'stdin')
    -o (--out)    output (can be 'stdout')
    -s (--step)   sliding window step (default: 5)
    -s (--size)   sliding window size (default: 60)

=head1 BUGS
  
=head1 REFERENCES
  
=head1 VERSION
  
  0.1
  
=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;
use POSIX qw/ceil floor/;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my ($size, $step) = (60, 5);
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "step|t=i" => \$step,
  "size|z=i" => \$size,
) or pod2usage(2);
pod2usage(1) if $help_flag;

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

my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');
my $cnt = 0;
while(my $seqO = $seqHI->next_seq()) {
  my ($id, $len, $seq) = ($seqO->id, $seqO->length, $seqO->seq);
  my ($beg, $end) = (1, $len);
  if($id =~ /^([\w\-]+)\:(\d+)\-(\d+)$/) {
    ($id, $beg, $end) = ($1, $2, $3);
  }
  
  my $n_win = int(($end-$beg+1-$size)/$step) + 1;
  next if $n_win < 1;

  my @wins;
  for my $i (0..$n_win-1) {
    my $beg1 = $beg + $step * $i;
    my $end1 = $beg + $step * $i + $size - 1;
    $end1 <= $end || die "$id:[$beg-$end] range error: [$beg1-$end1]\n";
    my $seq1 = substr($seq, $step*$i, $size);
    print $fho ">$id:$beg1-$end1\n";
    print $fho "$seq1\n";
  }
}
$seqHI->close();
close $fho;

exit 0;



