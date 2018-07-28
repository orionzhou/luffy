#!/usr/bin/env perl
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  fa2fq.pl - convert a Fasta file to Fastq file

=head1 SYNOPSIS
  
  fa2fq.pl [-help] [-in input] [-out output]

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
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my ($id, $seq, $len) = ("", "", 0);
while(<$fhi>) {
  chomp;
  if( /^\>(.+)/) {
    if($id ne "") {
      print $fho "\@$id\n";
      print $fho $seq."\n";
      print $fho "+\n";
      print $fho ("I" x $len) . "\n";
    }
    $id = $1;
    $seq = "";
    $len = 0;
  } else {
    $seq .= $_;
    $len += length($_);
  }
}
print $fho "\@$id\n";
print $fho $seq."\n";
print $fho "+\n";
print $fho ("I" x $len) . "\n";

close $fhi;
close $fho;



exit 0;
