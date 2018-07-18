#!/usr/bin/perl -w
use strict;
use DBI;
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy", {'RaiseError' => 1});
sub ideogram_localize()
{
  my ($chr,$bp) = @_;
  my $string_query = "select * from ideogram where chr='".$chr
    	."' and bp_start<".$bp." and bp_stop>".$bp ;
  my $sqr = $dbh->prepare($string_query);
  $sqr->execute();
  my $string = '';
  my $arm = '';
  while(my $ref = $sqr->fetchrow_hashref())
  {
    if(length($ref->{'band'}) > length($string))
    {
      $string = $ref->{'band'};
      $arm = $ref->{'arm'};
    }
  }
  return $chr.$arm.$string;
}