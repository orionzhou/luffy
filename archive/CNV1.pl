#!/usr/bin/perl -w
require 'chr_band.pl';
use strict;
use DBI;
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy", {'RaiseError' => 1});
my $file_prefix = "E:/Genius/papers/CNV/";
open( CNV, $file_prefix."CNV.txt" );

my @chr_arr = (1..22,"X","Y");
my $total = 0;
while(<CNV>)
{
  chomp;
  my @ele_arr = split("\t",$_);
  my $query  = "select * from DGV "
  	."where LocusID='".$ele_arr[1]."' and VariationType="
  	."'CopyNumber'";
  my $sqr = $dbh->prepare($query);
  $sqr->execute();
  my $count = 0;
  while(my $ref = $sqr->fetchrow_hashref())
  {
    print join("\t",substr($ref->{"Chr"},3,length($ref->{"Chr"})-1),
    	$ref->{"LocusID"},$ref->{"VariationID"},$ref->{"Landmark"},
    	$ref->{"Reference"},$ref->{"PubMedID"},$ref->{"Method"},$ref->{"Gain"},
        $ref->{"Loss"},$ref->{"TotalGainLossInv"},$ref->{"SampleSize"});
    print "\n";
    $count ++;
  }
  #print "$count\n";
}