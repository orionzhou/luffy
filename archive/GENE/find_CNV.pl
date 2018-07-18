#!/usr/bin/perl -w
require '../chr_band.pl';
use strict;
use DBI;
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy", {'RaiseError' => 1});
#SCZD my @gene_id_arr = ('84062','3084','267012','5625','1312','27185');
#Lipid my @gene_id_arr = ('336','116519','345','348','2194','3991','19');
#FMR1 
my @gene_id_arr = ('2332');
#DMD my @gene_id_arr = ('1756');
#DRD my @gene_id_arr = (50632,1812..1818);
foreach my $gene_id (@gene_id_arr)
{
  my $sqr = $dbh->prepare("SELECT * from gene where feature_id='GeneID:$gene_id' and ".
          "group_label='reference' and feature_type='GENE'");
  $sqr->execute();
  my $ref = $sqr->fetchrow_hashref();
  print join("\t",$ref->{'feature_name'},$ref->{'chr'},$ref->{'feature_id'},$ref->{'chr_start'},
        $ref->{'chr_stop'},$ref->{'chr_stop'}-$ref->{'chr_start'}+1,$ref->{'contig_orient'}
        ),"\n";
  my $chr = $ref->{'chr'};
  my $start = $ref->{'chr_start'};
  my $stop = $ref->{'chr_stop'};
  print "Up_1000kb\t+ gene +\tDown_1000k\n";
  #&find_CNV($chr,$start,$start);
  &find_CNV($chr,$start-1000000,$start+1000000);
  print "\n";
}
$dbh->disconnect();

sub find_CNV
{
  my ($chr,$start,$stop) = @_;
  my $query  = "select distinct(Locus_ID),Locus_Start,Locus_Stop from cnv "
  	."where Locus_Chr='chr$chr' and VariationType='CopyNumber'"
  	." and Locus_Stop>$start and Locus_Start<$stop";
  my $sqr = $dbh->prepare($query);
  $sqr->execute();
  my @locus_arr;
  while(my $ref = $sqr->fetchrow_hashref())
  {
    push(@locus_arr,join("\t",$ref->{'Locus_ID'},$ref->{'Locus_Start'},$ref->{'Locus_Stop'}
    	,ideogram_localize($chr,$start),ideogram_localize($chr,$stop)));
  }
  my $count = scalar(@locus_arr);
  print "$count";
  if($count>0)
  {
    foreach my $locus (@locus_arr)
    {
      print "\t$locus\n";
      my @temp = split("\t",$locus);
      #&fetch_detail_info($temp[0]);
    }
  }
}

sub fetch_detail_info
{
  my ($locus_ID) = @_;
  my $query  = "select * from cnv "
  	."where Locus_ID='$locus_ID' and VariationType='CopyNumber'";
  my $sqr = $dbh->prepare($query);
  $sqr->execute();
  while(my $ref = $sqr->fetchrow_hashref())
  {
    print "\t\t",join("\t",$ref->{'Variation_ID'},$ref->{'Landmark'},$ref->{'Start'}
    	,$ref->{'Stop'},$ref->{'PubMed_ID'},,$ref->{'Method'},$ref->{'Gain'}
        ,$ref->{'Loss'},$ref->{'Total'},$ref->{'Reference'},$ref->{'Sample_Size'});
  }
}