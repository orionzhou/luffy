#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 SYNOPSIS
  
  gtb.dedup.pl [-help] [-in input-file] [-out output-file]
    remove redundant gene models from a Gtb file (based on CDS locations)

  Options:
    -h (--help)   brief help message
    -i (--in)     input (Gtb) file
    -o (--out)    output (Gtb) file
    -f (--fam)    genefam priority file ($data/db/pfam/genefam.tbl)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Common;
use Location;
use Gtb;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fi, $fo) = ('') x 2;
my $ff = "$ENV{'data'}/db/pfam/genefam.tsv";
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "fam|f=s"  => \$ff,
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

my $hf;
open(my $fhf, "<$ff") or die "cannot read $ff\n";
while(<$fhf>) {
  chomp;
  /^fam/i && next;
  my ($fam, $dom, $pri) = split "\t";
  exists $hf->{$fam} && next;
  $hf->{$fam} = int($pri);
}
close $fhf;

my $t = readTable(-inh => $fhi, -header => 1);
close $fhi;

my $h;
for my $i (0..$t->lastRow) {
  my ($id, $par, $chr, $beg, $end, $srd, $locE, $locI, $locC, $loc5, $loc3, $phase, $src, $conf, $cat1, $cat2, $cat3, $note) = $t->row($i);
  my $rloc = [sort {$a->[0] <=> $b->[0]} @{locStr2Ary($locC)}];
  my $loc = $srd eq "-" ? [map {[$end-$_->[1]+1, $end-$_->[0]+1]} @$rloc] : 
    [map {[$beg+$_->[0]-1, $beg+$_->[1]-1]} @$rloc];
  my $str = join("|", $chr, locAry2Str($loc));
  my $len = $end - $beg + 1;
  $h->{$str} ||= [];
  push @{$h->{$str}}, $i;
}

my @idxsd;
for my $str (keys(%$h)) {
  my @idxs = @{$h->{$str}};
  @idxs > 1 || next;
 
  if(exists $hf->{$t->elm($idxs[0], "cat2")}) {
    @idxs = sort {$hf->{$t->elm($a, "cat2")} <=> $hf->{$t->elm($b, "cat2")}} @idxs;
  }
  my $idx = $idxs[0];
  push @idxsd, @idxs[1..$#idxs];
  
  my @notes = uniq( map {$t->elm($_, "note")} @idxs );
  $t->setElm($idx, "note", join(" ", @notes));
}
printf "%5d | %5d removed\n", scalar(@idxsd), $t->nofRow;

$t->delRows(\@idxsd);
print $fho $t->tsv(1);
close $fho;



__END__
