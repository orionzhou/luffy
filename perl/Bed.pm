package Bed;
use strict;
use Tabix;
use Data::Dumper;
use Common;
use Location;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/read_gap/;
@EXPORT_OK = qw//;

sub read_gap {
  my ($con, $chr, $beg, $end) = @_;
  my $iter = $con->query($chr, $beg - 1, $end);
  my @ary;
  return \@ary if ! $iter->get();
  while (my $line = $con->read($iter)) {
    my @ps = split("\t", $line);
    my ($chr, $beg, $end) = @ps;
    push @ary, [$chr, $beg, $end];
  }
  return \@ary;
}

1;
__END__
