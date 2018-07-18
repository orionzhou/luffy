#!/usr/bin/perl -w
use DBI;
use LWP::Simple;
my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
my $db     = "nucleotide";
my $report = "xml";

my $dbh = DBI->connect( "DBI:mysql:database=seq;host=localhost",
	"genius","prodigy", {'RaiseError' => 1});

#SCZD my @gene_id_arr = ('84062','3084','267012','5625','1312','27185');
#Lipid
 my @gene_id_arr = ('336','116519','345','348','2194','3991','19');
#FMR1 my @gene_id_arr = ('2332');
#DMD my @gene_id_arr = ('1756');

foreach $gene_id (@gene_id_arr)
{
  my $sqr = $dbh->prepare("SELECT * from gene where feature_id='GeneID:$gene_id' and ".
          "group_label='reference' and feature_type='GENE'");
  $sqr->execute();
  my $ref = $sqr->fetchrow_hashref();
  print join("\t",$ref->{'feature_name'},$ref->{'feature_id'},$ref->{'chr'},$ref->{'chr_start'},
  	$ref->{'chr_stop'},$ref->{'chr_stop'}-$ref->{'chr_start'}+1,$ref->{'contig_orient'}
  	),"\n";
  my $gene_start = $ref->{'chr_start'};
  my $gene_stop = $ref->{'chr_stop'};
  my $chr = $ref->{'chr'};
  #print "Upstream\n";
  #&find_str($chr, $gene_start-1000000, $gene_start-1);
  #print "InGene\n";
  #&find_str($chr, $gene_start, $gene_stop);
  #print "Downstream\n";
  #&find_str($chr, $gene_stop+1, $gene_stop+1000000);

  &find_str($chr, $gene_start-1000000, $gene_stop+1000000);
  #&find_str_strict($chr, $gene_start-1000000, $gene_stop+1000000);
  print "\n";
}
$dbh->disconnect();

sub find_str
{
    my ($chr,$str_start,$str_stop) = @_;
    my $query = "select distinct start from marshfield where chr_start>$str_start"
    	." and chr_stop<$str_stop and chr='$chr'";
    my $sqr1 = $dbh->prepare($query);
    $sqr1->execute();
    my $start_cM = 0;
    my $stop_cM = 0;
    while(my $ref0 = $sqr1->fetchrow_hashref())
    {
      my $ref_start = $ref0->{'start'};
      if($start_cM == 0 || $start_cM > $ref_start)
      {
        $start_cM = $ref_start;
      }
      if($stop_cM == 0 || $stop_cM < $ref_start)
      {
        $stop_cM = $ref_start;
      }
    }
    if($start_cM != 0)
    {
      $query = "select * from marshfield where start>$start_cM AND stop<$stop_cM AND chr='$chr'"
      	." OR start like '$start_cM%' OR start like '$stop_cM'";
      $sqr1 = $dbh->prepare($query);
      $sqr1->execute();
      while(my $ref1 = $sqr1->fetchrow_hashref())
      {
        print "\t",$ref1->{'name'};
      my $sqr = $dbh->prepare("select * from marker where name='".$ref1->{'name'}."'");
      $sqr->execute();
      my $genbanknum;
      if(my $ref = $sqr->fetchrow_hashref())
      {
        $genbanknum = $ref->{'genbanknum'};
        print "\t".$ref->{'alias'},"\t";
      }
      else
      {
        print "\t\t";
      }
      print join("\t",substr($ref1->{'sts_id'},0,length($ref1->{'sts_id'})-1),$ref1->{'chr_start'},
            $ref1->{'chr_stop'},$ref1->{'start'}),"\n";
      #&str_fetcher($genbanknum);
      }
    }
}

sub find_str_strict
{
    my ($chr,$str_start,$str_stop) = @_;
    my $query = "select * from marshfield where chr_start>$str_start"
    	." and chr_stop<$str_stop and chr='$chr'";
    my $sqr1 = $dbh->prepare($query);
    $sqr1->execute();
    while(my $ref1 = $sqr1->fetchrow_hashref())
    {
      print "\t",$ref1->{'name'};
      my $sqr = $dbh->prepare("select * from marker where name='".$ref1->{'name'}."'");
      $sqr->execute();
      my $genbanknum;
      if(my $ref = $sqr->fetchrow_hashref())
      {
        $genbanknum = $ref->{'genbanknum'};
        print "\t".$ref->{'alias'},"\t";
      }
      else
      {
        print "\t\t";
      }
      print join("\t",substr($ref1->{'sts_id'},0,length($ref1->{'sts_id'})-1),$ref1->{'chr_start'},
            $ref1->{'chr_stop'},$ref1->{'start'}),"\n";
    }
}

sub str_fetcher
{
  my ($genbanknum) = @_;
  if($genbanknum ne 'Unknown')
  {
    print "\t$genbanknum\t";
    my $query = "/efetch.fcgi?db=$db&id=$genbanknum&rettype=gp&retmode=$report";
    my $gf = get($utils.$query);
    while($gf =~ m|<GBKeyword>(.*)</GBKeyword>|g)
    {
      print $1."\t";
    }
    #$query = "/efetch.fcgi?db=$db&id=$genbanknum&rettype=gp&retmode=text";
    #$gf = get($utils.$query);
    #open(GB, ">E:/Genius/Gene/STS/$genbanknum.gb");
    #print GB $gf;
  }
  print "\n";
}