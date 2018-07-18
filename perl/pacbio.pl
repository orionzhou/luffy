#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  pacbio.pl - 

=head1 SYNOPSIS
  
  pacbio.pl [-help]

  Options:
    -h (--help)   brief help message
    -s (--stat)   print statistics

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Bio::SeqIO;
use Bio::AlignIO;
use Data::Dumper;
use Common;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $help_flag;
my $stat_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "stat|s"  => \$stat_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my @params = (
  ['HM034', 'HM034'],
  ['HM056', 'HM056'],
  ['HM340', 'HM340']
);
@params = (
  ['HM034', 'HM034.AC'],
  ['HM056', 'HM056.AC'],
  ['HM340', 'HM340.AC']
);

for (@params) {
  my ($qry, $tgt) = @$_;
  my $dir = "$ENV{'misc3'}/pacbio/$qry\_$tgt";
  -d $dir || make_path($dir);
  chdir $dir || die "cannot chdir to $dir\n";

  if($stat_flag) {
    exit;
  }
  my $f_tgt = "$ENV{'genome'}/$tgt/11_genome.fas";

  my $f_lst = "$dir/../reads/$qry.list";
  merge_seq_file($f_lst, "05.all.fastq");

  my $cmd = "blasr -nproc 16 -minReadLength 500 -minMatch 14 -bestn 1 \\
    -noSplitSubreads -clipping soft -sam \\
    -out 11.sam 05.all.fastq $f_tgt";
  runCmd($cmd);
  runCmd("samtools sort 11.sam 15");
  runCmd("samtools index 15.bam");
}

sub merge_seq_file {
  my ($fi, $fo) = @_;
  my @ary;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  while(<$fhi>) {
    chomp;
    /^\s*#/ && next;
    my $fp = $_;
    -s $fp || die "$fp not exist\n";
    push @ary, $fp;
  }
  if(@ary == 0) {
    print "no file in $fi\n";
  } elsif(@ary == 1) {
    runCmd("ln -sf $ary[0] $fo");
  } else {
    my $fi_str = join(" ", @ary);
    runCmd("cat $fi_str > $fo");
  }
}

__END__;
bedtools bamtofastq -i hm340_blasr/pb_blasr_ALPACA_sorted.bam -fq HM340_HM340/05.all.fastq
