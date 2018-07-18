#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  phy.denovo.pl - extract SNPs and construct phylogeny

=head1 SYNOPSIS
  
  phy.denovo.pl [-help]

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

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $dir = "$ENV{'misc3'}/phy_denovo";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my @orgs = qw/
  HM004 HM010 HM018 HM022 HM034
  HM050 HM056 HM058 HM060 HM095
  HM125 HM129 HM185 HM324 HM340
/;

#merge_snp_vcf(\@orgs, "00_raw", "01_allsnp.vcf");
#vcf_fill1(\@orgs, "01_allsnp.vcf", "02.bed", "03_split");
#vcf_fill2(\@orgs, "01_allsnp.vcf", "03_split", "05_filled.vcf");
#runCmd("bgzip -f 05_filled.vcf");
#runCmd("tabix -p vcf 05_filled.vcf.gz");

sub merge_snp_vcf {
  my ($orgs, $do, $fo) = @_;
  -d $do || make_path($do);

  for my $org (@$orgs) {
    my $fn = "$ENV{'misc3'}/$org\_HM101/23_blat/31.9/snp.vcf";
    runCmd("bgzip -cf $fn > $do/$org.vcf.gz");
    runCmd("tabix -p vcf $do/$org.vcf.gz");
  }
  my $fs_str = join(" ", map {"$do/$_.vcf.gz"} @orgs);
  runCmd("bcftools merge -o 01_allsnp.vcf $fs_str");
}
sub vcf_fill1 {
  my ($orgs, $fi, $fb, $do) = @_;
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

  -d $do || make_path($do);
  for my $org (@$orgs) {
    my $fg = "$ENV{'misc3'}/$org\_HM101/23_blat/31.9/gax.bed";
    runCmd("coverageBed -a $fg -b $fb | cut -f1,3,4 > $do/$org.bed");
  }
}
sub vcf_fill2 {
  my ($orgs, $fi, $do, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";

  my @fhbs;
  for my $org (@$orgs) {
    my $fb = "$do/$org.bed";
    open(my $fhb, "<$fb") or die "cannot read $fb\n";
    push @fhbs, $fhb;
  }
  while(<$fhi>) {
    chomp;
    if(/^\#/) {
      print $fho "$_\n";
      next;
    }
    my @ps = split "\t";
    @ps >= 10 || die "not 10 fileds: $_\n";
    my ($chr, $pos) = @ps[0,1];
    for my $idx (0..@$orgs-1) {
      my $gt = $ps[9+$idx];
      $gt =~ /^([\d\.])\/\1$/ or die "unknown gt: $gt\n";
      my $fhb = $fhbs[$idx];
      my $line = readline($fhb);
      
      next if $1 ne ".";
      chomp $line;
#      print $orgs[$idx]." ".$line."\n";
      my ($chr2, $pos2, $gt2) = split("\t", $line);
      ($chr eq $chr2 && $pos == $pos2) || 
        die "error $idx $chr:$pos $chr2:$pos2\n";
      $ps[9+$idx] = "0/0" if $gt2;
    }
    print $fho join("\t", @ps)."\n";
  }
  close $fhi;
  close $fho;
  for my $fhb (@fhbs) { close $fhb; };
}
my $d01 = "01_snp";
my $d11 = "11_snp_sample";
my $d12 = "12_phy";
my $d13 = "13_aln";
my $d21 = "21_phynj";
my $d22 = "22_phyml";
#make_path($d01, $d11, $d12, $d13, $d21, $d22);

my $reg = "chr5";
#runCmd("bcftools view -m2 -M2 -O u -v snps 05_filled.vcf.gz $reg | \\
#  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT]\\n' - | \\
#  snpfilter.pl -n 5 -m 1 -o 11.snp");

#runCmd("samplelines.pl -i 11.snp -o 12.snp -n 10000");
#runCmd("snp2phy.pl -i 12.snp -o 12.phy -r HM101");
#runCmd("phyml -i 12.phy -d nt");
#runCmd("mv 12.phy_phyml_tree.txt 12.nwk");
#runCmd("rm 12.phy_phyml*");

sub snp_convert {
    my ($dirI, $dirO, $chrs, $format) = @_;
    system("mkdir -p $dirO") unless -d $dirO;
    for my $chr (@$chrs) {
#    next if $chr eq "chr1";
        my $fo = "$dirO/$chr.$format";
        if($format =~ /^(aln)|(nexus)$/) {
            my $fi = "$dirI/$chr.phylip";
            runCmd("seqret $fi $fo -osformat $format");
        } elsif($format =~/^(phylip)|(structure)$/) {
            my $fi = "$dirI/$chr.ssp";
            runCmd("sspConvert -i $fi -o $fo -f $format");
        }
    }
}



