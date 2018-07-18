#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 SYNOPSIS
  
  gtb.pickalt.pl [-help] [-in input-file] [-out output-file]
      fix and pick alt model

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gtb)
    -o (--out)    output file (Gtb)

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
pod2usage(2) if !$fi || !$fo;

my ($fhi, $fho);
if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $t = readTable(-inh => $fhi, -header => 1);
close $fhi;

my $h;
my @idxsd;
for my $i (0..$t->lastRow) {
  my ($id, $par, $chr, $beg, $end, $srd, $locE, $locI, $locC, $loc5, $loc3, $phase, $src, $conf, $cat1, $cat2, $cat3, $note) = $t->row($i);
  if($cat1 ne "mRNA") {
    push @idxsd, $i;
    next;
  }
  my $loc = [sort {$a->[0] <=> $b->[0]} @{locStr2Ary($locC)}];
  my $len = locAryLen($loc);
  if(exists $h->{$par}) {
    my ($idxp, $lenp) = @{$h->{$par}};
    if($len <= $lenp) {
      push @idxsd, $i;
    } else {
      push @idxsd, $idxp;
      $h->{$par} = [$i, $len];
    }
  } else {
    $h->{$par} = [$i, $len];
  }
}
printf "%5d | %5d removed\n", scalar(@idxsd), $t->nofRow;

$t->delRows(\@idxsd);
print $fho $t->tsv(1);
close $fho;
printf "%5d models left\n", $t->nofRow;



__END__
