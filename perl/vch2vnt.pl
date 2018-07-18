#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  vch2vnt.pl - convert a VCF (human readable) file to SNP/IDM file

=head1 SYNOPSIS
  
  vch2vnt.pl [-help] [-in input-file] [-out output-prefix] 

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (VCH)
    -o (--out)    output prefix

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;
use Common;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;

my ($fi, $do) = ('') x 2;
my $fhi;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$do,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$do;

$fhi = \*STDIN;
unless ($fi eq "stdin" || $fi eq "-" || $fi eq "") {
  open ($fhi, "<$fi") || die "Can't open file $fi: $!\n";
}
-d $do || make_path($do);

my $hf = {
  'snp' => 1,
  'het' => 1, 
  'ins' => 1,
  'del' => 1,
};
for my $type (keys(%$hf)) {
  my $fn = "$do/$type"; 
  open(my $fh, ">$fn") || die "cannot write $fn\n";
  $hf->{$type} = $fh;
}

while( <$fhi> ) {
  chomp;
  my $line = $_;
  my @ps = split("\t", $line);
  @ps == 5 || die "not 5 lines:\n$line\n";
  my ($chr, $beg, $ref, $alts, $str) = @ps;
  my @alts = split(",", $alts);
  my $alt;
  if($str =~ /([A-Za-z0-9\-_]+)\=([\d\.])\/([\d\.])/) {
    if($2 eq $3) {
      $alt = $alts[$2-1];
    } elsif($2 ne $3 && $3 != 0) {
      my $alt1 = $2 == 0 ? $ref : $alts[$2-1];
      my $alt2 = $alts[$3-1];
      my $fh = $hf->{'het'};
      print $fh join("\t", $chr, $beg, $ref, "$alt1|$alt2")."\n";
      next;
    } else {
      die "unknown allele: $str\n";
    }
  } else {
    die "unknown allele: $str\n";
  }

  my ($lenr, $lena) = (length($ref), length($alt));
  my $fh;
  my ($tbeg, $tend, $talt);
  if($lenr == 1 && $lena == 1) {  
    $fh = $hf->{'snp'};
    print $fh join("\t", $chr, $beg, $ref, $alt)."\n";
  } elsif($lenr < $lena && $alt =~ /^$ref([ATCG]+)$/) {
    $fh = $hf->{'ins'};
    $talt = $1;
    $tbeg = $beg + $lenr - 1;
    print $fh join("\t", $chr, $tbeg, $tbeg, $talt)."\n";
  } elsif($lenr > $lena && $ref =~ /^$alt([ATCG]+)$/) {
    $fh = $hf->{'del'};
    $talt = $1;
    $tbeg = $beg + $lena;
    $tend = $beg + $lenr - 1;
    print $fh join("\t", $chr, $tbeg, $tend, $talt)."\n";
  } else {
    substr($ref, 0, 1) eq substr($alt, 0, 1) ||
      die "unknown allele: $ref $alt\n";
    $fh = $hf->{'del'};
    ($tbeg, $tend) = ($beg + 1, $beg + $lenr - 1);
    $talt = substr($ref, 1);
    print $fh join("\t", $chr, $tbeg, $tend, $talt)."\n";
    $fh = $hf->{'ins'};
    $tbeg = $beg;
    $talt = substr($alt, 1);
    print $fh join("\t", $chr, $tbeg, $tbeg, $talt)."\n";
  }
}
close $fhi;
for (keys(%$hf)) { close $hf->{$_}; }

runCmd("sort -k1,1 -k2,2n $do/snp -o $do/snp");
runCmd("bgzip -c $do/snp > $do/snp.gz");
runCmd("tabix -s 1 -b 2 -e 2 -f $do/snp.gz");

runCmd("sort -k1,1 -k2,2n -k3,3n $do/del -o $do/del");
runCmd("bgzip -c $do/del > $do/del.gz");
runCmd("tabix -s 1 -b 2 -e 3 -f $do/del.gz");

runCmd("sort -k1,1 -k2,2n $do/ins -o $do/ins");
runCmd("bgzip -c $do/ins > $do/ins.gz");
runCmd("tabix -s 1 -b 2 -e 3 -f $do/ins.gz");

