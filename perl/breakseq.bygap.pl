#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  breakseq.bygap.pl - breakup sequences by 'N's

=head1 SYNOPSIS
  
  breakseq.bygap.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (fasta)
    -o (--out)    output file (fasta)
    -g (--gap)    mininum gap length (def: 1000)

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
use FindBin;
use lib "$FindBin::Bin";
use Bio::SeqIO;
use Common;
use Location;

my ($fi, $fo) = ('') x 2;
my $min_gap = 1000;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "gap|g=i"  => \$min_gap,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my $base = [ split(/\./, $fo) ]->[0];
runCmd("fas2bed.pl -i $fi -o $base.1.bed");
runCmd("seqgap.pl -i $fi -o $base.2.gap.bed -m $min_gap");
runCmd("subtractBed -a $base.1.bed -b $base.2.gap.bed > $base.3.nogap.bed");
runCmd("seqret.pl -d $fi -b $base.3.nogap.bed -o $fo");

my @lens;
my $tb = readTable(-in => "$base.3.nogap.bed", -header => 0);
for my $i (0..$tb->lastRow) {
  my ($id, $beg, $end) = $tb->row($i);
  push @lens, $end - $beg;
}
@lens = reverse sort {$a <=> $b} @lens;
print "\n##### stat begins #####\n";
print join(" ", "max", @lens[0..4])."\n";
print join(" ", "min", reverse @lens[$#lens-4..$#lens])."\n";
print "##### stat ends   #####\n\n";
runCmd("rm $base.1.bed $base.2.gap.bed $base.3.nogap.bed");

exit 0;



