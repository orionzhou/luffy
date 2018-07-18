#!/usr/bin/perl -w
use DBI;
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=localhost",
	"genius","prodigy", {'RaiseError' => 1});
my $sqr = $dbh->prepare("select count(*) from marshfield");
$sqr->execute();
my $ref = $sqr->fetchrow_hashref();
my $total_count = $ref->{'count(*)'};
my $unfound_count = 0;
for(my $i=0; $i*200 <= $total_count; $i++)
{
  $sqr = $dbh->prepare("select * from marshfield where id>".($i*200)
  	." and id<".(($i+1)*200+1));
  $sqr->execute();
  while($ref = $sqr->fetchrow_hashref())
  {
    my $sts_id = substr($ref->{'sts_id'},0,length($ref->{'sts_id'})-1);
    my $sqr1 = $dbh->prepare("select chr_start,chr_stop from sts "
    	."where feature_id='$sts_id' and group_label='reference'");
    $sqr1->execute();
    if(my $ref1 = $sqr1->fetchrow_hashref())
    {
      my $string = "update marshfield set chr_start='".$ref1->{'chr_start'}
    	."',chr_stop='".$ref1->{'chr_stop'}."' where sts_id like '$sts_id%'";
      print "\tOkay\n";
      $sqr1 = $dbh->prepare($string);
      $sqr1->execute();
    }
    else
    {
      print ":(\n";
      $unfound_count ++;
    }
  }
}
print "\n$unfound_count not found";