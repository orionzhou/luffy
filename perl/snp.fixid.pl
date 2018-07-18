#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  snp.fixid.pl - fix NCGR non-standard HM*** ID

=head1 SYNOPSIS
  
  snp.fixid.pl [-help] [options] [-in input] [-out output] 

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my ($fhi, $fho);
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot read $fi\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $hm = {
  "HM017-I" => "HM017",
  "HM020-I" => "HM020",
  "HM022-I" => "HM022"
};

while( <$fhi> ) {
  chomp;
  my @ps = split "\t";
  my $nind = @ps - 4;
  for my $i (0..$nind-1) {
    if($ps[$i+4] =~ /([A-Za-z0-9\-_]+)\=([01\.])\/([01\.])/) {
      my ($oid, $nt1, $nt2) = ($1, $2, $3);
      my $id = exists $hm->{$oid} ? $hm->{$oid} : $oid;
      $ps[$i+4] = "$id=$nt1/$nt2";
    } else {
      die "unknown allele: $ps[$i+4]\n";
    }
  }
  print $fho join("\t", @ps)."\n";
}
close $fhi;
close $fho;


