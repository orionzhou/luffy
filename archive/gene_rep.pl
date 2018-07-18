#!/usr/bin/perl -w
use strict;
use DBI;
open(GENE,"F:/Genius/STR/X_STR/Xgene_final.txt");
open(REP,"F:/Genius/STR/X_STR/Xrepeat_final.txt");
open(WWW,">F:/Genius/STR/X_STR/gene_rep.txt");
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy");
my @gene_arr;
my @rep_arr;
while(<GENE>)
{
  chomp;
  @gene_arr = split("\t",$_);
  while(<REP>)
  {
    chomp;
    @rep_arr = split("\t",$_);
    if($rep_arr[1]<$gene_arr[1])
    {
      print WWW "\t",join("\t",@rep_arr,"\n");
    }
    else
    {
      print WWW join("\t","start",@gene_arr,"\n");
      my $length = length($_);
      seek(REP,-($length+2),1);
      last;
    }
  }
  while(<REP>)
  {
    chomp;
    @rep_arr = split("\t",$_);
    if($rep_arr[2]<$gene_arr[2])
    {
      print WWW "\t",join("\t",@rep_arr,"\n");
    }
    else
    {
      print WWW join("\t","stop",@gene_arr,"\n");
      my $length = length($_);
      seek(REP,-($length+2),1);
      last;
    }
  }
}