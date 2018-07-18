#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use File::Path qw/make_path remove_tree/;
use Eutils;
use Common;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my @orgs = qw/HM056 HM034 HM101 HM340/;

for my $org (@orgs) {
my $dir = "$ENV{'misc2'}/rnaseq/mt/31_cufflinks";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir $dir\n";

my $f_bam = "\$misc2/rnaseq/mt/22_tophat/$org\_$org/accepted_hits.bam";
my $f_gff = "\$genome/$org/51.gff";

make_path($org) unless -d $org;
my $cmd = sprintf("cufflinks -p %d -g $f_gff $f_bam -o $org", $ENV{'nproc'});
runCmd($cmd);
}

sub get_exp_support {
  my ($ff, $fi, $fo) = @_;
  
  my $tf = readTable(-in=>$ff, -header=>1);
  my $hf;
  for my $i (0..$tf->nofRow-1) {
    my ($id, $class, $idN, $idG, $name, $idT, $locus, $len, $cov, $fpkm, $fpkm_lo, $fpkm_hi, $fpkm_status) = $tf->row($i);
    $hf->{$id} = [$cov, $fpkm, $fpkm_status];
  }

  my $ti = readTable(-in=>$fi, -header=>1);
  open(FHO, ">$fo") or die "cannot write to $fo\n";
  print FHO join("\t", qw/id cov fpkm fpkm_status/)."\n";
  for my $i (0..$ti->nofRow-1) {
    my ($idQ, $tag, $idT) = $ti->row($i);
    my $stat = ["", "", ""];
    if($tag == 1) {
      if(exists $hf->{$idT}) {
        $stat = $hf->{$idT};
      } else {
        print "no exp for $idT\n";
      }
    } elsif($tag == 10) {
      next;
    } else {
      if(exists $hf->{$idQ}) {
        $stat = $hf->{$idQ};
      } else {
        print "no exp for $idQ\n";
      }
    }
    print FHO join("\t", $idQ, @$stat)."\n";
  }
  close FHO;
}


