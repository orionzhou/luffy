#!/usr/bin/perl -w
use DBI;
use LWP::Simple;
my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
my $db     = "SNP";
my $report = "FLT";

my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy", {'RaiseError' => 1});
my @chr_arr = (1..22,"X","Y");
my $file_prefix = "F:/Genius/SNP/";
foreach my $chr (@chr_arr)
{
  open( SNPs, ">".$file_prefix."Chr".$chr.".txt" ) or die("Couldn't create file!");
  &SNP_CHR($chr);
}
$dbh->disconnect();

sub SNP_CHR
{
  my ($chr) = @_;
  my $sqr = $dbh->prepare("SELECT * from gene where chromosome='$chr' and " .
	  "group_label='reference' and feature_type='GENE'");
  $sqr->execute();

  while(my $ref = $sqr->fetchrow_hashref())
  {
    print join("\t",$ref->{'id'},$ref->{'chr_start'},$ref->{'chr_stop'},$ref->{'chr_orient'},
    	$ref->{'contig_orient'},$ref->{'feature_name'},$ref->{'feature_id'}
        ),"\n";
    print SNPs join("\t",$ref->{'id'},$ref->{'chr_start'},$ref->{'chr_stop'},$ref->{'chr_orient'},
    	$ref->{'contig_orient'},$ref->{'feature_name'},$ref->{'feature_id'}
        ),"\n";
    my $snp_start;
    my $snp_stop;
    if($ref->{'contig_orient'} eq '+')
    {
      $snp_stop = $ref->{'chr_start'} - 1;
      $snp_start = $snp_stop - 5000 + 1;
    }
    elsif($ref->{'contig_orient'} eq '-')
    {
      $snp_start = $ref->{'chr_stop'} + 1;
      $snp_stop = $snp_start + 5000 - 1;
    }
    else
    {
      print "error!\n";
      exit;
    }
    my $query  = $ref->{'chromosome'} . "[CHR] AND $snp_start:$snp_stop"
    		. "[CHRPOS] AND human[ORGN]";
    my $esearch = "$utils/esearch.fcgi?" . "db=$db&retmax=1000&usehistory=y&term=";
    my $esearch_result = get($esearch . $query);
    $esearch_result =~ m|<Count>(\d+)</Count>.*<IdList>(.+)</IdList>|s;
    my $count = $1;
    my $ids = $2;
    if($count>0)
    {
      print "$count\n";
      print SNPs "$count\n";
      while( $ids =~ m|<Id>(\d+)</Id>|sg )
      {
        print $1,"\t";
        print SNPs $1,"\t";
      }
    }
    else
    {
      print "0";
      print SNPs "0";
    }
    print "\n\n";
    print SNPs "\n\n";
    #&esearch_result_handler($esearch_result);
  }
}

sub esearch_result_handler
{
  my ($esearch_result) = @_;
  $esearch_result =~
    m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;

  my $Count    = $1;
  my $QueryKey = $2;
  my $WebEnv   = $3;

  print "Count = $Count; QueryKey = $QueryKey; WebEnv = $WebEnv\n";

  my $retstart;
  my $retmax=3;

  for($retstart = 0; $retstart < $Count; $retstart += $retmax)
  {
    my $efetch = "$utils/efetch.fcgi?" .
                 "rettype=$report&retmode=text&retstart=$retstart&retmax=$retmax&" .
                 "db=$db&query_key=$QueryKey&WebEnv=$WebEnv";

    print "\nEF_QUERY=$efetch\n";

    my $efetch_result = get($efetch);

    print "---------\nEFETCH RESULT(".
           ($retstart + 1) . ".." . ($retstart + $retmax) . "): ".
          "[$efetch_result]\n";
  }
}