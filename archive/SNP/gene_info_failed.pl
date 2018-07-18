#!/usr/bin/perl -w
use strict;
use LWP::Simple;
my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
my $file_prefix = "E:/Genius/SNP/";
open (GENE, $file_prefix."SCZD_gene_2.txt" );
my @ele_arr;
while(<GENE>)
{
  chomp;
  @ele_arr = split("\t",$_);
  my $query  = $ele_arr[1]." AND human[ORGN]";
  my $esearch = "$utils/esearch.fcgi?db=UniGene&retmax=10&usehistory=y&term=";
  my $esearch_result = get($esearch . $query);
  $esearch_result =~ m|<Count>(\d+)</Count>.*<IdList>(.+)</IdList>|s;
  my $count = $1;
  my $ids = $2;
  if($count>0)
  {
    print "$count\n";
    while( $ids =~ m|<Id>(\d+)</Id>|sg )
    {
      print $1,"\t";
    }
  }
  else
  {
    print "0";
  }
  print "\n\n";
}