#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2gtf.pl - convert a Gtb file to GTF format

=head1 SYNOPSIS
  
  gtb2gtf.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file (Gtb)
      -out    output file (Gff)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Rna;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi;

my ($fhi, $fho);
if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $t = readTable(-in => $fi, -header => 1);
$t->sort("par", 1, 0, "id", 1, 0);
my $ref = group($t->colRef("par"));
my @gids = sort {$ref->{$a}->[0] <=> $ref->{$b}->[0]} keys %$ref;
printf "%d genes : %d models\n", scalar(@gids), $t->nofRow;

print $fho "##gff-version 3\n";
my ($cntG, $cntR) = (1, 1);
for my $gid (@gids) {
  my ($idxB, $cnt) = @{$ref->{$gid}};
  for my $i ($idxB..$idxB+$cnt-1) {
    my $rna = Rna->new(-gtb => $t->rowRef($i));
    print $fho $rna->to_gtf()."\n";
  }
}
close $fho;



__END__
