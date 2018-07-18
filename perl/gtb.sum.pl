#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.sum.pl - generate a simple report for a Gtb file

=head1 SYNOPSIS
  
  gtb.sum.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gtb)
    -o (--out)    output file (default: stdout)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gtb;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

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
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot read $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $t = readTable(-inh => $fhi, -header => 1);
close $fhi;

$t->sort("par", 1, 0, "id", 1, 0);
my @gids = $t->col("par");
my $ref = group(\@gids);
@gids = sort {$ref->{$a}->[0] <=> $ref->{$b}->[0]} keys %$ref;
printf $fho "%d genes : %d models\n", scalar(@gids), $t->nofRow;
my @cnts = uniq(map {$_->[1]} values %$ref);
for my $cnt (@cnts) {
  my @ids = grep {$ref->{$_}->[1] == $cnt} keys %$ref;
  printf $fho "  %6d x $cnt model(s)/gene [%s]\n", scalar(@ids), 
    @ids >= 0 ? $ids[0] : "";
}
close $fho;



__END__
