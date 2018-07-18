#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqdeltag.pl - calculate DeltaG free energe for a multi-fasta file

=head1 SYNOPSIS
  
  seqdeltag.pl [-help] [-in input-fasta] [-out output-tbl]

  Options:
      -help   brief help message
      -in     input fasta file (can be 'stdin')
      -out    output (can be 'stdout')

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Time::HiRes qw/gettimeofday tv_interval/;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

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

$ENV{"PATH"} = $ENV{"PATH"}.":/soft/unafold/3.8/bin";

my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');
my $t0 = [gettimeofday];
my $cnt = 0;

my $ftmp = "tmp";
while(my $seqO = $seqHI->next_seq()) {
    my ($id, $len, $seq) = ($seqO->id, $seqO->length, $seqO->seq);
    print $fho join("\t", $id, get_deltag($seq))."\n";
    
    $cnt ++;
    printf "%d: %.01f min\n", $cnt, tv_interval($t0, [gettimeofday]) / 60 if $cnt % 100000 == 0;
}
close $fhi;
close $fho;

sub get_deltag {
    my ($seq) = @_;
    my $dg;
    return $dg;
}

exit 0;
