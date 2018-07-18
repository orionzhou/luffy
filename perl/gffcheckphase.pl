#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gffphase.pl - check and fix the phases in a Gff file

=head1 SYNOPSIS
  
  gffphase.pl [-help] [-in input-Gff] [-seq refseq-fasta] [-out output-Gff]

  Options:
    -help   brief help message
    -in     input file (Gff)
    -out    output file (Gff)
    -seq    sequence-fasta

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gff;

my ($fi, $fo, $fs) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
    "seq|s=s" => \$fs,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fs;

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

print $fho "##gff-version 3\n";
my ($cntR, $cntG, $cntF) = (1, 1, 0);
my $it = parse_gff($fhi);
while(my $gene = $it->()) {
    print $fho $gene->to_gff()."\n";
    for my $rna ($gene->get_rna()) {
        $cntF += $rna->check_phase($fs);
        print $fho $rna->to_gff()."\n";
        printf "%5d RNA | %5d gene: %5d non-0 phase fixed\n", $cntR, $cntG, $cntF if $cntR % 1000 == 0;
        $cntR ++;
    }
    $cntG ++;
}

close $fhi;
close $fho;

__END__
