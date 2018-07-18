#!/usr/bin/perl -w
use DBI;
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy", {'RaiseError' => 1});

#SCZD my @gene_id_arr = ('84062','3084','267012','5625','1312','27185');
#Lipid my @gene_id_arr = ('336','116519','345','348','2194','3991','19');
#FMR1 my @gene_id_arr = ('2332');
#DMD
 my @gene_id_arr = ('1756','2332','215','959');

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
  print "Upstream\n";
  &find_ssr($chr, $gene_start-100000, $gene_start-1);
  print "InGene\n";
  &find_ssr($chr, $gene_start, $gene_stop);
  print "Downstream\n";
  &find_ssr($chr, $gene_stop+1, $gene_stop+100000);

  print "\n";
}
$dbh->disconnect();

sub find_ssr
{
    my ($chr,$str_start,$str_stop) = @_;
    my $query = "select * from SSR where chr_start>=$str_start"
    	." and chr_stop<=$str_stop and chr='$chr' and feature_name REGEXP '^\\\\([A-Z]{3,4}\\\\)n\$'";
    my $sqr = $dbh->prepare($query);
    $sqr->execute();
    while(my $ref = $sqr->fetchrow_hashref())
    {
      my $ssr_start = $ref->{'chr_start'};
      my $ssr_stop = $ref->{'chr_stop'};
      my $ssr_id = $ref->{'id'};
      my $ssr_link = "http://202.117.161.53/app/contents/show_ssr.php?id=".$ssr_id;
      print "\t";
      print join("\t",$ssr_id,$ref->{'feature_name'},$ssr_start,$ssr_stop,$ssr_link);
      print "\n";
    }
}