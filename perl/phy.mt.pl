#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  phy.mt.pl - extract SNPs and construct phylogeny

=head1 SYNOPSIS
  
  phy.mt.pl [-help] [-opt option]

  Options:
    -h (--help)   brief help message
    -t (--opt)    option

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

my $opt = 'pan16x';
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "opt|t=s" => \$opt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$opt;

my $dir = "$ENV{'misc1'}/phy.mt/$opt";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";
my $fv = "../acc319.vcf.gz";
-s $fv || die "$fv not there\n";

my $sams = get_sample_ids("01.ids.txt");
my $n_sam = scalar(@$sams);
my $sam = join(",", @$sams);
print "==========================\n";
printf "opt[%s] %d samples:\n%s\n", $opt, $n_sam, $sam;

my ($nap_max, $mac_min) = (0.2, 1);
my $na_max = int($n_sam * $nap_max);
print "==========================\n";
printf "missing data <= %d\n", $na_max;
printf "minor allele count >= %d\n", $mac_min;
print "==========================\n";

my $reg = ""; #"chr5";
runCmd("bcftools view -U -m2 -M2 -O u -v snps -s $sam $fv $reg | \\
  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT]\\n' - | \\
  snpfilter.pl -n $na_max -m $mac_min | \\
  snp.fixid.pl -o 11.snp");

runCmd("samplelines.pl -i 11.snp -o 31.snp -n 100000");
runCmd("snp2phy.pl -i 31.snp -o 31.phy");

runCmd("\$soft/phyml/phyml -i 31.phy -d nt");
runCmd("mv 31.phy_phyml_tree.txt 31.nwk");
runCmd("rm 31.phy_phyml*");

sub get_sample_ids {
  my ($fp) = @_;
  open(my $fhi, "<$fp") or die "cannot read $fp\n";
  
  my $hm = {
    "HM017" => "HM017-I",
    "HM020" => "HM020-I",
    "HM022" => "HM022-I"
  };
  
  my @sams;
  while(<$fhi>) {
    chomp;
    while(/(HM[\d\-\_I]+)/ig) {
      my $sam = $1;
      $sam = $hm->{$sam} if exists $hm->{$sam};
      push @sams, $sam;
    }
  }
  return \@sams;
}

##runCmd("seqret $f12 $f13 -osformat aln");
##runCmd("clustalw2 -INFILE=$f13 -BOOTSTRAP=1000 -OUTORDER=INPUT -OUTPUTTREE=phylip -BOOTLABELS=node -CLUSTERING=NJ -KIMURA");
##runCmd("mv $d13/$reg.phb $d21");

__END__
sub snp_convert {
  my ($dirI, $dirO, $chrs, $format) = @_;
  system("mkdir -p $dirO") unless -d $dirO;
  for my $chr (@$chrs) {
#  next if $chr eq "chr1";
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

