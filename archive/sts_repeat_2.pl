#!/usr/bin/perl -w
use DBI;
use strict;
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius", "prodigy", {'RaiseError'=>1});
open(RRR,"F:/Genius/STR/X_STR/sts_DXS_rep.txt");
my @sts_arr;
my @sts_repeat;
my $repeat_count = 0;
my $sts_count = 0;
while(<RRR>)
{
  chomp;
  @sts_arr = split("\t",$_);
  if(scalar(@sts_arr)>5)
  {
    my $i = 0;
    while($sts_arr[5+3*$i])
    {
      print join("\t",$sts_arr[0],@sts_arr[5+3*$i..7+3*$i],"\n");
      my $sql = $dbh->prepare("select repeat_id from Xsts where id=".$sts_arr[0]);
      $sql->execute();
      my $ref = $sql->fetchrow_hashref();
      my $repeat_id = $ref->{'repeat_id'};
      if($repeat_id ne '')
      {
        $repeat_id .= "*".$sts_arr[5+3*$i];
      }
      else
      {
        $repeat_id = $sts_arr[5+3*$i];
      }
      $sql = $dbh->prepare("update Xsts set repeat_id='".$repeat_id
      	."' where id=".$sts_arr[0]);
      $sql->execute();
      print "\t$repeat_id";

      $sql = $dbh->prepare("select sts_id from Xrepeat where id=".$sts_arr[5+3*$i]);
      $sql->execute();
      $ref = $sql->fetchrow_hashref();
      my $sts_id = $ref->{'sts_id'};
      if($sts_id ne '')
      {
        $sts_id .= "*".$sts_arr[0];
      }
      else
      {
        $sts_id = $sts_arr[0];
      }
      $sql = $dbh->prepare("update Xrepeat set sts_id='".$sts_id
      	."' where id=".$sts_arr[5+3*$i]);
      $sql->execute();
      print "\t$sts_id\n";

      $i++;
      $repeat_count++;
    }
  $sts_count++;
  }
}
print "$sts_count STS were mapped to $repeat_count repeats";