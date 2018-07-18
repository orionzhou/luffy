#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  repeatmasker.pl - run RepeatMasker and parse output

=head1 SYNOPSIS
  
  repeatmasker.pl [-help] [-org genome-id]

  Options:
    -h (--help)   brief help message
    -g (--org)    genome ID

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;

my ($org) = ('');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "org|g=s"  => \$org,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$org;

my $dir = "$ENV{'genome'}/$org";
chdir $dir or die "cannot chdir to $dir\n";
my $fg = "11_genome.fas";

-d "12.rm" && runCmd("rm -rf 12.rm");
mkdir "12.rm";
runCmd("RepeatMasker -pa 24 -species medicago -xsmall -dir $dir/12.rm $dir/$fg");
parse_rm("12.rm/11_genome.fas.out", "12.rm.tbl");
runCmd("awk 'BEGIN{OFS=\"\\t\"} {print \$1, \$2-1, \$3, \$9 \" | \" (\$5)}' 12.rm.tbl > 12.rm.raw.bed");
runCmd("sortBed -i 12.rm.raw.bed > 12.rm.bed");
runCmd("bedToBigBed -tab -type=bed4 12.rm.bed $fg 12.rm.bb");
runCmd("rm 12.rm.raw.bed");

sub parse_rm {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  while(<$fhi>) {
    chomp;
    my @ps = split " ";
    next if @ps == 0 || $ps[0] !~ /^[0-9e\.]+/i;
    my ($score, $div, $del, $ins, $qid, $qbeg, $qend, $qleft,
      $srd, $tid, $tfam, $tbeg, $tend, $tleft, $id) = @ps;
    if($srd eq "C") {
      $srd = "-";
      ($tbeg, $tleft) = ($tleft, $tbeg);
    }
    print $fho join("\t", $qid, $qbeg, $qend, $id, 
      $tid, $tbeg, $tend, $srd, $tfam, $score, $div, $del, $ins)."\n";
  }
  close $fhi;
  close $fho;
}


exit 0;
