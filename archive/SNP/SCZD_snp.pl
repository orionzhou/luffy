#!/usr/bin/perl -w
use strict;
use DBI;
my $file_prefix = "E:/Genius/SNP/";
open (GENE, $file_prefix."SCZD_gene_3.txt" );

my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy", {'RaiseError' => 1});
my @ele_arr;
while(<GENE>)
{
  chomp;
  @ele_arr = split("\t",$_);
  print join("\t",@ele_arr[1,6],"\t");
  my $query  = "select * from gene_snp where GeneID='".$ele_arr[6]."'";
  my $sqr = $dbh->prepare($query);
  $sqr->execute();
  if(my $ref = $sqr->fetchrow_hashref())
  {
    print join("\t",$ref->{'snp_count'},$ref->{'snp'}),"\n";
  }
  else
  {
    print "\n";
  }
}