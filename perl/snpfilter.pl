#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  snpfilter.pl - read, filter and write Snp file

=head1 SYNOPSIS
  
  snpfilter.pl [-help] [options] [-in input] [-out output] 

  Options:
    -h, --help       brief help message
    -i, --in         input
    -o, --out        output
    -n, --missing    samples of missing data no more than ?
    -m, --minor-ac   minor allele count no less than ?

=head1 VERSION
  
  0.1
  
=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my ($co_mis, $co_mac) = (-1, 0);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
    "missing|n=i" => \$co_mis,
    "minor-ac|m=i" => \$co_mac,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my ($fhi, $fho);
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot read $fi\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my @names;
my ($nind, $npos) = (0, 0);

my $first = 1;
while( <$fhi> ) {
  chomp;
  my $line = $_;
  my @ps = split("\t", $line);
  if($first) {
    $nind = @ps - 4;
  } else {
    die "not $nind + 4 cols\n" unless $nind + 4 == @ps;
  }
  my ($n_anc, $n_der, $n_mis) = (0, 0, 0);
  for my $i (0..$nind-1) {
    if($ps[$i+4] =~ /([A-Za-z0-9\-_]+)\=([01\.])\/([01\.])/) {
      push @names, $1 if $first;
      if($2 eq "0" && $3 eq "0") {
        $n_anc ++;
      } elsif($2 eq "1" && $3 eq "1") {
        $n_der ++;
      } else {
        $n_mis ++;
      }
    } else {
      die "unknown allele: $ps[$i+4]\n";
    }
  }
  my $mac = min($n_anc, $n_der);
  
  my $pass = 1;
  $pass = 0 if $co_mis > -1 && $n_mis > $co_mis;
  $pass = 0 if $mac < $co_mac;

  print $fho $line."\n" if $pass;
  $first = 0 if $first;
}
close $fhi;
close $fho;


sub sspExpand {
  my ($dirI, $dirO, $ids, $chrs) = @_;
  for my $chr (@$chrs) {
#    next unless $chr eq "chr1";
    my $fi = file($dirI, "$chr.txt");
    for my $id (@$ids) {
      system("mkdir -p $dirO/$id") unless -d "$dirO/$id";
    }
    runCmd("sspConvert -i $fi -o $dirO -f plain");
  }
}
sub ssp2Tbl {
  my ($dirI, $dirO, $ids, $chrs, $refDb) = @_;
  system("mkdir -p $dirO") unless -d $dirO;
  my @fns;
  for my $chr (@$chrs) {
    my $fi = file($dirI, "$chr.txt");
    my $fo = file($dirO, "$chr.tbl");
    runCmd("sspConvert -i $fi -o $fo -f table");
    runCmd('sed -i \'s/^\([0-9\.][0-9\.]*\)/'.$chr.'\t\1\t\1/\' '.$fo);
    push @fns, $fo;
  }
  my $fo = file($dirO, "01.tbl.gz");
  runCmd("cat ".join(" ", @fns)." | bgzip > $fo");
  runCmd("tabix -f -s1 -b2 -e3 $fo");
  runCmd("rm ".join(" ", @fns));
}
sub ssp2Ped {
  my ($dirI, $dirO, $chrs, $refDb) = @_;
  system("mkdir -p $dirO") unless -d $dirO;
  for my $chr (@$chrs) {
    my $fi = file($dirI, "$chr.txt");
    my $fo = file($dirO, "$chr.ped");
    runCmd("sspFilter -i $fi -o $fo -b 1 -f ped");
  }
}
sub snpEff {
  my ($d01, $d31, $d32, $chrs) = @_;
  system("mkdir -p $d31") unless -d $d31;
  system("mkdir -p $d32") unless -d $d32;
#  $chrs = [map { "chr".$_ } (4)];
  for my $chr (@$chrs) {
    my $f01 = "$d01/$chr.txt";
    my $f31 = "$d31/$chr.vcf";
    my $f32 = "$d32/$chr.tbl";
#    runCmd("ssp2Vcf -i $f01 -o $f31 -r HM000 -c $chr");
    my $cmd = "java -Xmx5g -jar \$src/snpEff_2_0_3/snpEff.jar eff mt_35 -c \$m/conf/snpEff.config -i vcf -ud 0 $f31 -s $d32/$chr.html > $f32";
    system($cmd);
  }
}


