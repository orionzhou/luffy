#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  phy.finderror.pl - 

=head1 SYNOPSIS
  
  phy.finderror.pl [-help]

  Options:
    -h (--help)   brief help message

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
use List::MoreUtils qw/first_index/;

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $dir = "$ENV{'misc3'}/phy_finderror";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $fv = "$ENV{'misc3'}/hapmap/30_vnt/acc319.vcf.gz";
my $fc = "$ENV{'misc3'}/hapmap/30_vnt/ingroup.txt";
#runCmd("bcftools view -U -m2 -M2 -c2 -O z -v snps -S $fc -o 01.vcf.gz $fv chr5");

my @orgs = ("HM018B", "HM018C", "HM017B", "HM022B");
for my $org (@orgs) {
  my $dirc = "$ENV{'misc3'}/$org\_HM101/23_blat/31.9";
#  runCmd("cp -f $dirc/snp.vcf 07.$org.vcf");
}
# make manual corrections

for my $org (@orgs) {
  my $dirc = "$ENV{'misc3'}/$org\_HM101/23_blat/31.9";
#  runCmd("bgzip -cf 07.$org.vcf > 07.$org.vcf.gz");
#  runCmd("tabix -f -p vcf 07.$org.vcf.gz"); 
#  runCmd("bcftools view -O z -o 08.$org.vcf.gz 07.$org.vcf.gz chr5");
#  runCmd("tabix -f -p vcf 08.$org.vcf.gz");
}

my $vcfstr = join(" ", map {"08.$_.vcf.gz"} @orgs);
#runCmd("bcftools merge 01.vcf.gz $vcfstr -o 11.vcf");
#vcf2bed("11.vcf", "11.bed");
for my $org (@orgs) {
  my $dirc = "$ENV{'misc3'}/$org\_HM101/23_blat/31.9";
#  runCmd("coverageBed -a $dirc/gax.bed -b 11.bed | cut -f1,3,4 > 15.$org.bed");
}

my $fcs = [ map {"15.$_.bed"} @orgs ];
#vcf_fill("11.vcf", \@orgs, $fcs, "21.vcf");
#runCmd("bgzip -f 21.vcf");
#runCmd("tabix -f -p vcf 21.vcf.gz");

sub vcf2bed {
  my ($fi, $fb) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fhb, ">$fb") or die "cannot write $fb\n";
  while(<$fhi>) {
    chomp;
    /^\#/ && next;
    my @ps = split "\t";
    @ps >= 10 || die "not 10 fields: $_\n";
    my ($chr, $pos) = @ps[0..1];
    print $fhb join("\t", $chr, $pos-1, $pos)."\n";
  }
  close $fhi;
  close $fhb;
}
sub vcf_fill {
  my ($fi, $orgs, $fcs, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";

  my @fhcs;
  for my $fc (@$fcs) {
    open(my $fhc, "<$fc") or die "cannot read $fc\n";
    push @fhcs, $fhc;
  }
  my @cidxs;
  while(<$fhi>) {
    chomp;
    if(/^\#/) {
      if(/^\#CHROM/) {
        my @cnames = split("\t", $_);
        for my $org (@$orgs) {
          my $cidx = first_index {$_ eq $org} @cnames;
          $cidx >= 0 || die "cannot find $org in header\n";
          push @cidxs, $cidx;
        }
      }
      print $fho "$_\n";
      next;
    }
    my @ps = split "\t";
    @ps >= 10 || die "not 10 fileds: $_\n";
    my ($chr, $pos) = @ps[0,1];
    
    for my $i (0..$#cidxs) {
      my $gt = $ps[$cidxs[$i]];
      $gt =~ /^([\d\.])\/\1/ or die "unknown gt: $gt\n";
      my $al = $1;
    
      my $line = readline($fhcs[$i]);
      chomp $line;
#    print $orgs[$idx]." ".$line."\n";
      my ($chr2, $pos2, $cov) = split("\t", $line);
      die "error $i $chr:$pos $chr2:$pos2\n" unless ($chr eq $chr2 && $pos == $pos2);
      if($al eq "." && $cov) {
        $ps[$cidxs[$i]] = "0/0";
      }
    }
    print $fho join("\t", @ps)."\n";
  }
  close $fhi;
  close $fho;
}

#runCmd("bcftools view -m2 -M2 -O u -v snps 21.vcf.gz chr5 | \\
#  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT]\\n' - | \\
#  snpfilter.pl -n 80 -m 1 -o 51.snp");

runCmd("samplelines.pl -i 51.snp -o 52.snp -n 11000");
runCmd("snp2phy.pl -i 52.snp -o 52.phy -r HM101_ref");
runCmd("phyml -i 52.phy -d nt");
runCmd("mv 52.phy_phyml_tree.txt 52.nwk");
runCmd("rm 52.phy_phyml*");




