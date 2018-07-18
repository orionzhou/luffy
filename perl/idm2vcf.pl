#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  idm2vcf.pl - convert small InDels to VCF file

=head1 SYNOPSIS
  
  idm2vcf.pl [-help] [-qry query-genome] [-tgt target-genome]

  Options:
    -h (--help)   brief help message
    -q (--qry)    query genome (def: HM056)
    -t (--tgt)    target genome (def: HM101)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Common;
use Location;
use Seq;
use Vcfhead;
use List::Util qw/min max sum/;

my ($qry, $tgt) = ('HM056', 'HM101');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "qry|q=s" => \$qry,
  "tgt|t=s" => \$tgt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$qry || !$tgt;

my $data = $ENV{'data'};
my $qry_fas = "$data/genome/$qry/11_genome.fas";
my $tgt_fas = "$data/genome/$tgt/11_genome.fas";
my $qdb = Bio::DB::Fasta->new($qry_fas);
my $tdb = Bio::DB::Fasta->new($tgt_fas);

my $dir = "$data/misc3/$qry\_$tgt/23_blat";

my $fi = "$dir/31.9/idm";
my $fo = "$dir/31.9/idm.vcf";

open(my $fhi, "<$fi") or die "cannot read $fi\n";
open(my $fho, ">$fo") or die "cannot write $fo\n";

print $fho $vcfhead."\n";
print $fho "#".join("\t", @colhead, $qry)."\n";

while( <$fhi> ) {
  chomp;
  next if /(^\#)|(^id)|(^\s*$)/s;
  my ($tid, $tbeg, $tend, $tsrd, $qid, $qbeg, $qend, $qsrd, $cid, $lev) 
    = split "\t";
  $lev == 1 || next;
  my $tlen = $tend - $tbeg - 1;
  my $qlen = $qend - $qbeg - 1;
  ($tlen < 50 && $qlen < 50) || next;
  $tsrd eq "+" or die "$tid:$tbeg-$tend not +\n";
  my $ref = seqret_simple($tdb, $tid, $tbeg, $tend - 1, "+");
  my $alt = $qsrd eq "+" ? 
    seqret_simple($qdb, $qid, $qbeg, $qend - 1, $qsrd) :
    seqret_simple($qdb, $qid, $qbeg + 1, $qend, $qsrd);
  if($ref eq $alt) {
#    print "$tid:$tbeg-$tend $ref $qid:$qbeg-$qend $alt\n";
    next;
  }
#  substr($ref, 0, 1) eq substr($alt, 0, 1) || print "$tid:$tbeg-$tend ref[$ref] alt[$alt]\n";
  print $fho join("\t", $tid, $tbeg, ".", $ref, $alt, 50, '.',
    '.', 'GT', '1/1')."\n";
}
close $fhi;
close $fho;


__END__
