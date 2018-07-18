#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  mt.nbs.pl - 

=head1 SYNOPSIS
  
  mt.nbs.pl [-help] [-org organism]

  Options:
    -h (--help)   brief help message
    -g (--org)    genmome ID (organism) to process

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Location;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($org) = ('');
my $help_flag;
#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "org|g=s"  => \$org,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$org;

my $dir = "$ENV{'genome'}/$org/42.nbs";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir $dir\n";

my $fg = "../augustus/31.gtb";
-s $fg || die "$fg is not there\n";

my $fm = "$ENV{'misc2'}/nbs/mt_40/30.all.hmm";

runCmd("hmmsearch.pl -i ../augustus/31.fas -m $fm -o 01.tbl");
build_nbs_gtb("01.tbl", $fg, "11.gtb");
runCmd("gtb.pickalt.pl -i 11.gtb -o 12.gtb");
if($org eq "HM1010") {
  runCmd("liftover.gene.pl -i 12.gtb -o 16.liftover.tbl \\
    -r $ENV{'genome'}/HM101/11_genome.fas \\
    -g $ENV{'genome'}/HM101/51.gtb \\
    -x $ENV{'misc3'}/$org\_HM101/23_blat/41.9/gax.gz \\
    -s $ENV{'misc3'}/$org\_HM101/23_blat/41.9/snp.gz \\
    -d $ENV{'misc3'}/$org\_HM101/23_blat/41.9/idm.gz");
}

sub build_nbs_gtb {
  my ($fi, $fg, $fo) = @_;
  my $ti = readTable(-in => $fi, -header => 1);
  my $tg = readTable(-in => $fg, -header => 1);
  my %h = map {$ti->elm($_, "hid") => $ti->elm($_, "qid")} 0..$ti->lastRow;
  
  my @idxs;
  for my $i (0..$tg->lastRow) {
    my ($id, $par, $chr, $beg, $end, $srd, $locE, $locI, $locC, $loc5, $loc3, $phase, $src, $conf, $cat1, $cat2, $cat3, $note) = $tg->row($i);
    exists $h{$id} || next;
    push @idxs, $i;
    $cat3 = $h{$id};
    $cat2 = '';
    if($cat3 =~ /^tnl/i) {
      $cat2 = "TIR-NBS-LRR";
    } else {
      $cat3 =~ /^cnl/i || die "unknown cat: $cat3\n";
      $cat2 = "CC-NBS-LRR";
    }
    $tg->setElm($i, "cat2", $cat2);
    $tg->setElm($i, "cat3", $cat3);
    $h{$id} = 0;
  }
  printf "found %d | %d ids\n", scalar(@idxs), $ti->nofRow;
  my @ids_missed = grep {$h{$_} != 0} keys(%h);
  if(@ids_missed > 0) {
    print "  IDs not found: ".join(" ", @ids_missed)."\n";
  }
  
  my $to = $tg->subTable(\@idxs, [$tg->header]);
  open(my $fho, ">$fo") or die "cannot read $fo\n";
  print $fho $to->tsv(1);
  close $fho;
}


