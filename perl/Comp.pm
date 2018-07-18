package Comp;
use strict;
use Tabix;
use Data::Dumper;
use Common;
use Location;
use Seq;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/read_mcl partition_record/;
@EXPORT_OK = qw//;


sub read_mcl {
  my ($fi) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  my @groups;
  while(<$fhi>) {
    chomp;
    my @ids = split "\t";
    push @groups, \@ids;
  }
  close $fhi;
  return \@groups;
}
sub partition_record {
  my ($groups, $orgs) = @_;
  my @ary;
  for (@$groups) {
    my $h = { map {$_ => []} @$orgs };
    my $hn = { map {$_ => 0} @$orgs };
    my @ids = @$_;
    for (@ids) {
      my ($org, $gid) = split /\-/;
      push @{$h->{$org}}, $gid;
      $hn->{$org} ++;
    }

    my $n = max(values(%$hn));
    for my $i (1..$n) {
      my @ort;
      for my $org (@$orgs) {
        if($hn->{$org} >= $i) {
          push @ort, $h->{$org}->[$i-1];
        } else {
          push @ort, '';
        }
      }
      push @ary, \@ort;
    }
  }
  return \@ary;
}


1;
__END__
