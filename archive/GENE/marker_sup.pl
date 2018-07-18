#!/usr/bin/perl -w
use DBI;
use LWP::Simple;
my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
my $db     = "nucleotide";
my $report = "xml";

my $dbh = DBI->connect( "DBI:mysql:database=seq;host=localhost",
	"genius","prodigy", {'RaiseError' => 1});
my $sqr = $dbh->prepare("select * from marker where genbanknum='Unknown'");
$sqr->execute();
my $update_count = 0;
my @id_remained;
while(my $ref = $sqr->fetchrow_hashref())
{
  print join("\t",$ref->{'name'},$ref->{'alias'},"\n");
  my $esearch_query = "/esearch.fcgi?db=$db&retmax=20&term=".$ref->{'name'};
  if($ref->{'alias'} ne "Unknown")
  {
    $esearch_query .= " OR ".$ref->{'alias'};
  }
  my $esearch_result = get($utils.$esearch_query);
  $esearch_result =~ m|<Count>(\d+)</Count>.*<IdList>(.*)</IdList>|s;
  my $count = $1;
  my $ids = $2;
  my @id_arr;
  print $count,"\t";
  if($count > 0)
  {
    while($ids =~ m|<Id>(\d+)</Id>|g)
    {
      print $1,"\t";
      push(@id_arr,$1);
    }
    print "\n";
    my $esum_query = "/esummary.fcgi?db=$db&id=".join(",",@id_arr)
       ."&rettype=gp&retmode=$report";
    my $esum_result = get($utils.$esum_query);
    $esum_result =~ m|<Id>(\d+)</Id>(\s)<Item Name="Caption" Type="String">([^N]\d+)</Item>|;
    if($3)
    {
      print "\t$1\t$3\n";
      my $sqr1 = $dbh->prepare("update marker set genbanknum='$3' where name='"
      	.$ref->{'name'}."'");
      $sqr1->execute();
      $update_count ++;
    }
    else
    {
      push(@id_remained,$ref->{'id'});
    }
  }
  else
  {
    print "\n";
  }
  print "\n";
}
print "$update_count items updated:)\n\n";
print join("\n",@id_remained);