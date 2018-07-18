#!/usr/bin/perl -w
use strict;
my @chr_arr = (1..22,"X","Y");
my $file_prefix = "E:/Genius/SNP/SNP_by_Chr/";
my $total_gene_count = 0;
my $total_SNP_count = 0;
foreach my $chr (@chr_arr)
{
  open( SNPs, $file_prefix."Chr".$chr.".txt" ) or die("Couldn't open file");
  my $temp;
  my $p2;
  my $gene_count = 0;
  my $snp_count = 0;
  while(!eof(SNPs))
  {
    $temp = <SNPs>;
    $p2 = <SNPs>;
    if($p2 ne "\n")
    {
      $temp = <SNPs>;
      $snp_count += int($p2);
    }
    $temp = <SNPs>;
    $gene_count ++;
  }
  print "Chromosome $chr\t$gene_count\t$snp_count\n";
  $total_gene_count += $gene_count;
  $total_SNP_count += $snp_count;
}
print "\t$total_gene_count\t$total_SNP_count";