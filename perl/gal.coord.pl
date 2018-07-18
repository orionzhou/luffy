#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal.coord.pl - transform coordinate system

=head1 SYNOPSIS
  
  gal.coord.pl [-help] [-in input-file] [-opt option] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -o (--opt)    option ('qry' or 'tgt', default: 'qry')
    -q (--qry)    qry-size file [required if opt='qry']
    -t (--tgt)    tgt-size file [required if opt='tgt']

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gal;

my ($fi, $fo) = ('') x 2;
my ($fq, $ft) = ('') x 2;
my $opt = "qry";
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "opt|p=s"  => \$opt,
  "qry|q=s"  => \$fq,
  "tgt|t=s"  => \$ft
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$opt;
$opt = lc($opt);

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

my ($hq, $ht);
if($opt eq "qry") {
  $hq = read_size($fq);
} elsif($opt eq "tgt") {
  $hq = read_size($ft);
} else {
  die "unknown opt: $opt\n";
}

print $fho join("\t", @HEAD_GAL)."\n";
while(<$fhi>) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my @ps = split "\t";
  next unless @ps == 21;
  if($opt eq 'qry') {
    my ($qId, $qBeg, $qEnd) = @ps[6..8];
    if($qId =~ /^(\w+)\-([0-9e\+]+)\-([0-9e\+]+)$/) {
      my ($qi, $beg, $end) = ($1, $2, $3);
      $ps[6] = $qi;
      $ps[7] = $beg + $qBeg - 1;
      $ps[8] = $beg + $qEnd - 1;
      exists $hq->{$qi} || die "no size for $qi\n";
      $ps[10] = $hq->{$qi};
    }
  } elsif($opt eq 'tgt') {
    my ($tId, $tBeg, $tEnd) = @ps[1..3];
    if($tId =~ /^(\w+)\-([0-9e\+]+)\-([0-9e\+]+)$/) {
      my ($ti, $beg, $end) = ($1, $2, $3);
      $ps[1] = $ti;
      $ps[2] = $beg + $tBeg - 1;
      $ps[3] = $beg + $tEnd - 1;
      exists $ht->{$ti} || die "no size for $ti\n";
      $ps[5] = $ht->{$ti};
    }
  }
  print $fho join("\t", @ps)."\n";
}
close $fhi;
close $fho;

sub read_size {
  my ($fs) = @_;
  open(my $fhs, "<$fs") or die "cannot read $fs\n";
  my $h;
  while(<$fhs>) {
    chomp;
    next if /^\s*$/;
    my ($id, $size) = split "\t";
    !exists $h->{$id} || die "$id seen twice in $fs\n";
    $h->{$id} = $size;
  }
  return $h;
}


__END__
