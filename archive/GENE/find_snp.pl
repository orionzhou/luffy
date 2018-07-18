#!/usr/bin/perl -w
require "chr_band.pl";
use DBI;
use LWP::Simple;
my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
my $db     = "SNP";
my $report = "FLT";

my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy", {'RaiseError' => 1});
#SCZD my @gene_id_arr = ('84062','3084','267012','5625','1312','27185');
#Lipid my @gene_id_arr = ('336','116519','345','348','2194','3991','19');
#FMR1 my @gene_id_arr = ('2332');
#DMD my @gene_id_arr = ('1756');
#HTR2A my @gene_id_arr = ('3356');
#DRD my @gene_id_arr = (50632,1812..1818);
#my @gene_id_arr = (51639,3575,977,283685,116444,3150,10084);
#DRD1-5 my @gene_id_arr = (1812..1816);
#OPRM1 my @gene_id_arr = (4988);
#INSR
 my @gene_id_arr = (3643);

foreach $gene_id (@gene_id_arr)
{
  #&find_gene_info($gene_id);
  &find_snp($gene_id,1,10000);
  &find_snp($gene_id,2);
  &find_snp($gene_id,3,5000);
}
$dbh->disconnect();

sub find_snp
{
  my ($gene_id,$choice,$length) = @_;
  my $sqr = $dbh->prepare("SELECT * from gene where feature_id='GeneID:$gene_id' and ".
	  "group_label='reference' and feature_type='GENE'");
  $sqr->execute();
  my $ref = $sqr->fetchrow_hashref();
  print join("\t",$ref->{'feature_name'},$ref->{'feature_id'},$ref->{'chr'},$ref->{'chr_start'},
    	$ref->{'chr_stop'},$ref->{'chr_stop'}-$ref->{'chr_start'}+1,$ref->{'contig_orient'}
        ),"\n";

  my $snp_start;
  my $snp_stop;
  if($ref->{'contig_orient'} eq '+')
  {
      if($choice == 1)
      {
        $snp_stop = $ref->{'chr_start'} - 1;
        $snp_start = $snp_stop - $length + 1;
      }
      elsif($choice == 2)
      {
        $snp_start = $ref->{'chr_start'};
        $snp_stop = $ref->{'chr_stop'};
      }
      elsif($choice == 3)
      {
        $snp_start = $ref->{'chr_stop'} + 1;
        $snp_stop = $snp_start + $length - 1;
      }
      else
      {
        print "error!\n";
        exit;
      }
  }
  elsif($ref->{'contig_orient'} eq '-')
  {
      if($choice == 1)
      {
        $snp_start = $ref->{'chr_stop'} + 1;
        $snp_stop = $snp_start + $length - 1;
      }
      elsif($choice == 2)
      {
        $snp_start = $ref->{'chr_start'};
        $snp_stop = $ref->{'chr_stop'};
      }
      elsif($choice == 3)
      {
        $snp_stop = $ref->{'chr_start'} - 1;
        $snp_start = $snp_stop - $length + 1;
      }
      else
      {
        print "error!\n";
        exit;
      }
  }
  else
  {
      print "error!\n";
      exit;
  }
  if($choice == 1)
  {
      print "Upstream ".$length."bp\t";
  }
  elsif($choice == 2)
  {
      print "Inside of Gene\t";
  }
  elsif($choice == 3)
  {
      print "Downstream ".$length."bp\t";
  }
  &fetch_snp_brief($ref->{'chr'},$snp_start,$snp_stop);
  #&fetch_snp_detail($ref->{'chr'},$snp_start,$snp_stop);
  print "\n";
}

sub fetch_snp_brief
{

  my ($chr,$snp_start,$snp_stop) = @_;
  my $query  = $chr . "[CHR] AND $snp_start:$snp_stop"
    		. "[CHRPOS] AND human[ORGN]";
  my $esearch = "$utils/esearch.fcgi?" . "db=$db&retmax=100&usehistory=y&term=";
  my $esearch_result = get($esearch . $query);
  $esearch_result =~
    m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;

  my $Count    = $1;
  my $QueryKey = $2;
  my $WebEnv   = $3;

  if($Count>0)
  {
      print "$Count\n";
      #while( $ids =~ m|<Id>(\d+)</Id>|sg )
      #{
      #  print $1,"\t";
      #}
    #print "Count = $Count; QueryKey = $QueryKey; WebEnv = $WebEnv\n";
    my $retstart;
    my $retmax=20;
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
  else
  {
      print "0\n";
  }
}

sub fetch_snp_detail
{
  my ($chr,$snp_start,$snp_stop) = @_;
  my $query  = $chr . "[CHR] AND $snp_start:$snp_stop"
    		. "[CHRPOS] AND human[ORGN]";
  my $esearch = "$utils/esearch.fcgi?" . "db=$db&retmax=100&usehistory=y&term=";
  my $esearch_result = get($esearch . $query);
  $esearch_result =~
    m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;

  my $Count    = $1;
  my $QueryKey = $2;
  my $WebEnv   = $3;

  if($Count>0)
  {
      print "$Count\n";
      #while( $ids =~ m|<Id>(\d+)</Id>|sg )
      #{
      #  print $1,"\t";
      #}
    #print "Count = $Count; QueryKey = $QueryKey; WebEnv = $WebEnv\n";
    my $retstart;
    my $retmax=20;
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
        print "\t$result\n";
      }
    }
  }
  else
  {
      print "0\n";
  }
}

sub find_gene_info
{
  my ($gene_id) = @_;
  my $sqr = $dbh->prepare("SELECT * from gene where feature_id='GeneID:$gene_id' and ".
	  "group_label='reference'");
  $sqr->execute();
  my $ref = $sqr->fetchrow_hashref();
  my $start_band = ideogram_localize($ref->{'chr'},$ref->{'chr_start'});
  my $stop_band = ideogram_localize($ref->{'chr'},,$ref->{'chr_stop'});
  print join("\t",$ref->{'feature_name'},$ref->{'feature_id'},$ref->{'feature_type'},
        $ref->{'chr'},
        $ref->{'chr_start'},$ref->{'chr_stop'},$ref->{'chr_stop'}-$ref->{'chr_start'}+1,
        $ref->{'contig_orient'},$start_band,$stop_band),"\n";
  while($ref = $sqr->fetchrow_hashref())
  {
    print join("\t",$ref->{'feature_name'},$ref->{'feature_id'},$ref->{'feature_type'},
    	$ref->{'chr'},
    	$ref->{'chr_start'},$ref->{'chr_stop'},$ref->{'chr_stop'}-$ref->{'chr_start'}+1,
        $ref->{'contig_orient'}),"\n";
  }
  print "\n";
}