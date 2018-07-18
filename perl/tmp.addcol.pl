#!/usr/bin/perl -w

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Location;
use Gal;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('', '');
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $head_flag = 1;
while( <$fhi> ) {
  chomp;
  my @ps = split "\t";
  @ps == 20 || die "not 20 fields\n";
  my $col = "";
  if($head_flag) {
    $col = "lev";
    $head_flag = 0;
  } 
  print $fho join("\t", @ps[0..10], $col, @ps[11..19])."\n";
}
close $fhi;
close $fho;

__END__
