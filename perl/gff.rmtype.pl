#!/usr/bin/perl -w
=pod BEGIN
  
=head1 NAME
  
  gff.rmtypes.pl: remove features of certein type(s) with their children from a GFF3 file

=head1 SYNOPSIS
  
 gffrmtypes.pl [-help] -t <type> -i <input-file> -o <output-file>

 Options:
   -h (--help)    brief help message
   -i (--in)      input file
   -o (--out)     output file
   -t (--type)    types to remove

=cut

use strict;
use FindBin;
use lib "$FindBin::Bin";

use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw/make_path/;
use Common;
use Gff;

my $help_flag;
my ($fi, $fo) = ('', '');
my $types = '';
my %options = (
  "help|h" => \$help_flag,
  "in|i=s" => \$fi,
  "out|o=s" => \$fo,
  "type|t=s" => \$types,
);

#--------------------------------- MAIN -----------------------------------#
GetOptions(%options) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$types;
my %ht = map {$_ => 1} split(",", $types);

my ($fhi, $fho);
if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my %hi;
while(<$fhi>) {
  chomp;
  if(/^\#/) {
    print $fho $_."\n";
  } else {
    my @ps = split "\t";
    my ($chr, $src, $type, $beg, $end, $score, $srd, $phase, $tagstr) = @ps;
    my $htag = parse_gff_tags($tagstr);
    
    my $skip = 0;
    if(exists $ht{$type}) {
      $hi{$htag->{'ID'}} = 1;
      $skip = 1;
    }
    if(exists $htag->{'Parent'} && exists $hi{$htag->{'Parent'}}) {
      $hi{$htag->{'ID'}} = 1 if exists $htag->{'ID'};
      $skip = 1;
    }
    print $fho join("\t", @ps)."\n" if $skip == 0;
  }
}
close $fhi;
close $fho;


