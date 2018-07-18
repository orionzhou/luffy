#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  pfam2tbl.pl - convert Pfam(TSV) file to Tbl file

=head1 SYNOPSIS
  
  pfam2tbl.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (pfam tsv) file
    -o (--out)    output (tbl) file
    -e            maximum E-value (def: 1)
    -l (--len)    minimum length (def: 10)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Location;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fi, $fo) = ('') x 2;
my ($max_e, $min_len) = (1, 10);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "e=f"     => \$max_e,
  "len|l=i" => \$min_len,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi;

my ($fhi, $fho);
if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my ($h, $hn);
while(<$fhi>) {
  chomp;
  my ($id, $md5, $len, $src, $fam, $note, $beg, $end, $e) = split "\t";
  $e <= $max_e || next;
  $h->{$id} ||= [];
  push @{$h->{$id}}, [$beg, $end, $fam, $e];
  $hn->{$fam} ||= $note;
  $hn->{$fam} eq $note || die "$fam has >1 notes\n";
}
close $fhi;

print $fho join("\t", qw/id beg end fam e note/)."\n";
for my $id (sort(keys(%$h))) {
  my @ary = @{$h->{$id}};
  my $locs = [ map {[$_->[0], $_->[1]]} @ary ];
  my $stats = [ map {$_->[3]} @ary ];
  my $ref = tiling($locs, $stats, 1);

  my (@fams, @begs, @ends, $es);
  for (@$ref) {
    my ($beg, $end, $idx) = @$_;
    $end - $beg + 1 >= $min_len || next;
    my ($fam, $e) = ($ary[$idx]->[2], $ary[$idx]->[3]);
    print $fho join("\t", $id, $beg, $end, $fam, $e, $hn->{$fam})."\n";
  }
}
close $fho;

exit 0;
