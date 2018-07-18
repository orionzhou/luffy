#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  rename.pl - rename fields with a map file

=head1 SYNOPSIS
  
  rename.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -m (--map)    ID mapping file
    -c (--col)    columns to rename (e.g., '6,12'; default: '1')
    -s (--skip)   skip lines with less than # cols (default: 1)

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
use List::Util qw/min max sum/;

my ($fi, $fo, $fm) = ('', '', '');
my $colstr = "1";
my $colskip = 1;
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "map|m=s"  => \$fm,
  "col|c=s"  => \$colstr,
  "skip|s=i" => \$colskip,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fm;

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

my $h;
open(my $fhm, "<$fm") or die "cannot read $fm\n";
while(<$fhm>) {
  chomp;
  my @ps = split "\t";
  my ($oid, $nid) = @ps[0, -1];
  exists $h->{$oid} && die "$oid has >1 mappings\n";
  $h->{$oid} = $nid;
}
close $fhm;

my @cols = split(",", $colstr);
while( <$fhi> ) {
  chomp;
  my $line = $_;
  if($line =~ /(^id)|(^\#)|(^\s*$)/) {
    print $fho $line."\n";
    next;
  }
  my @ps = split("\t", $line, -1);
  if(scalar(@ps) < $colskip) {
    print $fho join("\t", @ps)."\n";
    next;
  }
  for my $col (@cols) {
    my $oid = $ps[$col - 1];
    exists $h->{$oid} || die "no mapping for $oid\n";
    $ps[$col - 1] = $h->{$oid};
  }
  print $fho join("\t", @ps)."\n";
}
close $fhi;
close $fho;

__END__
