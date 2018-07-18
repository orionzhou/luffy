#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use Net::SFTP;

my $ftp = Net::SFTP->new("data.ncgr.org", Debug=>0) || die "failed connect\n";
$ftp->login("pipestone", "smRRkn++") || die "login failed ", $ftp->message;
$ftp->cwd("vcf_Mt4.0") || die $ftp->message;

my $dir = "/home/youngn/zhoup/Data/misc3/hapmap_mt40/30_vnt";

my @fns = qw/
  all_lines.unified_genotyper_round2.bqsr.filter_dp_recode_gt.vcf.gz/;
for my $fn (@fns) {
  $ftp->get($fn, "$dir/$fn") || die $ftp->message;
}
