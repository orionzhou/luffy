#!/usr/bin/perl -w
use strict;
my @chr_arr = (1..22,"X","Y");
my $file_prefix = "E:/Genius/SNP/SNP_by_Chr/";
open( FILE, ">".$file_prefix."SNPs.txt" );
foreach my $chr (@chr_arr)
{
  open( SNPs, $file_prefix."Chr".$chr.".txt" ) or die("Couldn't open file");
  my $p1;
  my $p2;
  my $p3;
  my $temp;
  my @ele_arr_1;
  my @ele_arr_3;
  my $snp_count;
  while(!eof(SNPs))
  {
    $p1 = <SNPs>;
    chomp($p1);
    @ele_arr_1 = split("\t",$p1);
    $snp_count = 0;
    $p2 = <SNPs>;
    $p3 = "";
    if($p2 ne "\n")
    {
      $p3 = <SNPs>;
      $snp_count = int($p2);
    }
    $temp = <SNPs>;
    chomp($p3);
    @ele_arr_3 = split("\t",$p3);
    $temp = join("#",@ele_arr_3);
    print FILE join("\t",@ele_arr_1[0,5,6],$snp_count,$temp),"\n";
  }
}