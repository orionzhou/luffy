#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  msf2aln.pl - convert alignment file format (MSF -> ALN)

=head1 SYNOPSIS
  
  msf2aln.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input alignment (MSF)
    -o (--out)    output alignment file (ALN)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::AlignIO;

#--------------------------------- MAIN -----------------------------------#
my ($fi, $fo) = ('') x 2;
my $help_flag;
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my $alni = Bio::AlignIO->new(-file => "<$fi", -format => "msf");
my $alno = Bio::AlignIO->new(-file => ">$fo", -format => "clustalw");
while(my $aln = $alni->next_aln) {
  $alno->write_aln($aln);
}

exit 0;
