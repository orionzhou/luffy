#!/usr/bin/perl -w
use LWP::Simple;
my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
my $db     = "SNP";
my $report = "FLT";

#COMT my @location = ('22','18320069','18330069');
#NRG1 my @location = ('8','32614282','32624282');
#DMD
#my @location = ('X','33139594','33149594');
#my @location = ('X','33056466','33066466');
#my @location = ('X','32948238','32958238');
#my @location = ('X','32340292','32350292');
#my @location = ('X','32083507','32093507');
#my @location = ('X','31436336','31446336');
my @location = ('X','31194945','31204945');

&find_snp_upstream_gene(@location);

sub find_snp_upstream_gene
{
  my ($chr,$start,$stop) = @_;
  my $query  = $chr."[CHR] AND $start:$stop"."[CHRPOS] AND human[ORGN]";
  my $esearch = "$utils/esearch.fcgi?" . "db=$db&retmax=1000&usehistory=y&term=";
  my $esearch_result = get($esearch . $query);
  $esearch_result =~ m|<Count>(\d+)</Count>.*<IdList>(.+)</IdList>|s;
  my $count = $1;
  my $ids = $2;
  if($count>0)
  {
      print "$count\n";
      #while( $ids =~ m|<Id>(\d+)</Id>|sg )
      #{
      #  print $1,"\t";
      #}
  }
  else
  {
    print "0";
  }
  #&fetch_snp_detail($esearch_result);
  &fetch_snp_brief($esearch_result);
  print "\n";
}

sub fetch_snp_detail
{
  my ($esearch_result) = @_;
  $esearch_result =~
    m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;

  my $Count    = $1;
  my $QueryKey = $2;
  my $WebEnv   = $3;
  #print "Count = $Count; QueryKey = $QueryKey; WebEnv = $WebEnv\n";
  my $retstart;
  my $retmax=3;
  for($retstart = 0; $retstart < $Count; $retstart += $retmax)
  {
    my $efetch = "$utils/efetch.fcgi?" .
                 "rettype=$report&retmode=text&retstart=$retstart&retmax=$retmax&" .
                 "db=$db&query_key=$QueryKey&WebEnv=$WebEnv";
    my $efetch_result = get($efetch);
    my @result_arr = split("\n\n",$efetch_result);
    foreach $result (@result_arr)
    {
      $result =~ s|^\n||;
      $result =~ s|\s.\s|\t|g;
      $result =~ s|\n|\n\t|g;
      print "\t$result\n";
    }
    print "\n";
  }
}

sub fetch_snp_brief
{
  my ($esearch_result) = @_;
  $esearch_result =~
    m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;

  my $Count    = $1;
  my $QueryKey = $2;
  my $WebEnv   = $3;
  #print "Count = $Count; QueryKey = $QueryKey; WebEnv = $WebEnv\n";
  my $retstart;
  my $retmax=3;
  for($retstart = 0; $retstart < $Count; $retstart += $retmax)
  {
    my $efetch = "$utils/efetch.fcgi?" .
                 "rettype=$report&retmode=text&retstart=$retstart&retmax=$retmax&" .
                 "db=$db&query_key=$QueryKey&WebEnv=$WebEnv";
    my $efetch_result = get($efetch);
    my @result_arr = split("\n\n\n\n",$efetch_result);
    foreach $result (@result_arr)
    {
      $result =~ s|^\n||;
      $result =~ s|\s.\s|\t|g;
      $result =~ s|\n|\n\t|g;
      $result =~ m|\trs(\d+)|;
      print "\trs$1\t";
      $result =~ m|CTG.*reference.*chr\-pos=(\d+)|;
      print "$1\n";
    }
  }
}