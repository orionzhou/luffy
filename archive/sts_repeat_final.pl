#!/usr/bin/perl -w
use DBI;
use strict;
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius", "prodigy", {'RaiseError'=>1});
open(RRR,"F:/Genius/STR/X_STR/sts_DXS.txt");
open(WWW,">F:/Genius/STR/X_STR/sts_DXS_all.txt");
my @sts_arr;
my @sts_repeat;
my $repeat_count = 0;
my $sts_count = 0;
while(<RRR>)
{
  chomp;
  @sts_arr = split("\t",$_);
  print WWW join("\t",@sts_arr[0..3]);
  my $sql = $dbh->prepare("select * from Xrepeat_all where chr_start >= "
  	.$sts_arr[1]." and chr_stop <= ".$sts_arr[2]);
  $sql->execute();
  while( my $ref = $sql->fetchrow_hashref() )
  {
  	print WWW "\t",join("\t",$ref->{'id'},$ref->{'chr_start'},$ref->{'chr_stop'},
        	$ref->{'feature_name'});
        $repeat_count++;
  }
  print "\n";
  $sts_count++;
}
print "$sts_count STS were mapped to $repeat_count repeats";