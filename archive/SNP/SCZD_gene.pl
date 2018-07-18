#!/usr/bin/perl -w
use strict;
my $file_prefix = "E:/Genius/SNP/";
open( GENE, $file_prefix."SCZD_gene.txt" );
my @ele_arr;
while(<GENE>)
{
  chomp;
  if($_ =~ /^\*(\d+)\s(.*);\s([a-zA-Z0-9]+)$/ )
  {
    print join("\t",$1,$3,$2),"\n";
  }
  elsif($_ =~ /^\*(\d+)\s(.*)$/)
  {
    my @name_arr = split(" ",$2);
    print join("\t",$1,$2,$name_arr[0]),"\n";
  }
}