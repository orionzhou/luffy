#!/usr/bin/perl -w
use strict;
use DBI;
my $file_prefix = "E:/Genius/SNP/";
open (GENE, $file_prefix."SCZD_gene_2.txt" );

my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy", {'RaiseError' => 1});
my @ele_arr;
while(<GENE>)
{
  chomp;
  @ele_arr = split("\t",$_);
  print $ele_arr[1],"\t";
  my $query  = "select * from gene where feature_name='".$ele_arr[1]
  	."' and feature_type='GENE' and group_label='reference'";
  my $sqr = $dbh->prepare($query);
  $sqr->execute();
  if(my $ref = $sqr->fetchrow_hashref())
  {
    print join("\t",$ref->{'chromosome'},$ref->{'chr_start'},$ref->{'chr_stop'},
    	$ref->{'feature_id'});
    my $ideogram = &ideogram_localize($ref->{'chromosome'},$ref->{'chr_start'});
    print "\t",$ideogram,"\n";
  }
  else
  {
    print "\n";
  }
}

sub ideogram_localize()
{
  my ($chr,$bp) = @_;
  my $string_query = "select * from ideogram where chr='".$chr
    	."' and bp_start<".$bp." and bp_stop>".$bp
        ." and band regexp '^[0-9]{2}\$'" ;
  my $sqr = $dbh->prepare($string_query);
  $sqr->execute();
  my $ref = $sqr->fetchrow_hashref();
  return $ref->{'chr'}.$ref->{'arm'}.$ref->{'band'};
}