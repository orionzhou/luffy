#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  blast2psl.pl - convert a BLAST tabular output to PSL format

=head1 SYNOPSIS
  
  blast2psl.pl [-help] [-in input-file] [-out output-file]
    -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen length nident mismatch gaps evalue bitscore qseq sseq'

  Options:
      -help   brief help message
      -in     input file
      -out    output file

=head1 DESCRIPTION

  This program converts an input BLAST tabular file to a PSL file 

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Blast;

my ($fi, $fo) = ('') x 2;
my ($fhi, $fho);
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
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
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "stdout" || $fo eq "-") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

while(<$fhi>) {
    chomp;
    my $ps = [ split "\t" ];
    next unless @$ps == 16;
    $ps = blast2psl($ps);
    print $fho join("\t", @$ps)."\n";
}
__END__
