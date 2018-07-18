#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  parse.rm.pl

=head1 SYNOPSIS
  
  parse.rm.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file 
    -o (--out)    output file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;

my ($fi, $fo) = ('') x 2;
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

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
while(<$fhi>) {
  chomp;
  my @ps = split " ";
  next if @ps == 0 || $ps[0] !~ /^[0-9e\.]+/i;
  my ($score, $div, $del, $ins, $qid, $qbeg, $qend, $qleft,
    $srd, $tid, $tfam, $tbeg, $tend, $tleft, $id) = @ps;
  if($srd eq "C") {
    $srd = "-";
    ($tbeg, $tleft) = ($tleft, $tbeg);
  }
  print $fho join("\t", $qid, $qbeg, $qend, $id, 
    $tid, $tbeg, $tend, $srd, $tfam, $score, $div, $del, $ins)."\n";
}
close $fhi;
close $fho;

__END__
