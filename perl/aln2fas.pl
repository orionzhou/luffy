#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  aln2fas.pl - convert Aln (clustalw) format to Fasta format

=head1 SYNOPSIS
  
  aln2fas.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (aln) file
    -o (--out)    output (fasta) file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::AlignIO;

my ($fi, $fo) = ('') x 2;
my $len = 10;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "len|l=i" => \$len,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my $ai = Bio::AlignIO->new( -file => "<$fi", -format => "clustalw" );
my $ao = Bio::AlignIO->new( -file => ">$fo", -format => "fasta" );

while( my $aln = $ai->next_aln() ) {
  $ao->write_aln($aln);
}

exit 0;
