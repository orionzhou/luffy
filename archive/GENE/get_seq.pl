#!/usr/bin/perl -w
require "chr_band.pl";
use DBI;
use LWP::Simple;
my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
my $db     = "nucleotide";
my $retmode = "text";

my $dbh = DBI->connect( "DBI:mysql:database=seq;host=localhost",
	"genius","prodigy", {'RaiseError' => 1});
my $file_prefix = "E:/Genius/Gene/sequence/";

#SCZD
 my @gene_id_arr = ('84062','3084','267012','5625','1312','27185');
#Lipid my @gene_id_arr = ('336','116519','345','348','2194','3991','19');
#FMR1 my @gene_id_arr = ('2332');
#DMD my @gene_id_arr = ('1756');
foreach $gene_id (@gene_id_arr)
{
  my $sqr = $dbh->prepare("SELECT * from gene where feature_id='GeneID:$gene_id' and ".
	  "group_label='reference' and feature_type='GENE'");
  $sqr->execute();
  my $ref = $sqr->fetchrow_hashref();
  my $start;
  my $stop;
  if($ref->{'contig_orient'} eq '+')
  {
      $stop = $ref->{'contig_start'} - 1;
      $start = $stop - 10000 + 1;
  }
  elsif($ref->{'contig_orient'} eq '-')
  {
      $start = $ref->{'contig_stop'} + 1;
      $stop = $start + 10000 - 1;
  }
  else
  {
      print "error!\n";
      exit;
  }
  &get_seq($ref->{'contig'},$start,$stop,$ref->{'feature_name'});
}

sub get_seq
{
  my ($id,$start,$stop,$gene) = @_;
  my $query = "/efetch.fcgi?db=$db&id=$id&seq_start=$start&seq_stop=$stop"
  	."&rettype=fasta&retmode=$retmode";
  open(FAS, ">".$file_prefix.$gene."_up_10k.fas");
  my $buffer = get($utils.$query);
  print $buffer;
  print FAS $buffer;
}