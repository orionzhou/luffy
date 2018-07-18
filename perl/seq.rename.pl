#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seq.rename.pl - order & rename seqs by descending seq-length

=head1 SYNOPSIS
  
  seq.rename.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (must be a file)
    -o (--out)    output file (default: stdout)
    -p (--pre)    ID prefix (default: 'scf')

=cut
  
#### END of POD documentation.
#------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Common;

my ($fi, $fo) = ('') x 2;
my $pre = "scf";
my $help_flag;

#------------------------------ MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "pre|p=s" => \$pre,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi;

my $fho;
if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

runCmd("seqlen.py $fi $fi.sizes");

my $h;
open(my $fhs, "<$fi.sizes") or die "cannot read $fi.sizes\n";
while(<$fhs>) {
  chomp;
  my ($id, $len) = split "\t";
  ! exists $h->{$id} || die "2 seqs with same ID: $id\n";
  $h->{$id} = $len;
}
close $fhs;

my $db = Bio::DB::Fasta->new($fi);
my @ids = sort {$h->{$b} <=> $h->{$a} || $a cmp $b} keys(%$h);
my $dig = getDigits(scalar(@ids));

my $seqHO = Bio::SeqIO->new(-fh => $fho, -format => 'fasta');
open(my $fhm, ">$fi.map") or die "cannot write $fi.map\n";
for my $i (0..$#ids) {
  my $id = $ids[$i];
  my $nid = sprintf "$pre%0".$dig."d", $i;
  
  my $seq = $db->seq($id);
  
  $seqHO->write_seq( Bio::Seq->new(-id => $nid, -seq => $seq) );
  print $fhm join("\t", $id, $h->{$id}, $nid)."\n";
}
close $fhm;
$seqHO->close();

runCmd("rm $fi.sizes $fi.index");

exit 0;



