#!/usr/bin/perl -w
use strict;
use DBI;
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy", {'RaiseError' => 1});

# query
my $sqr = $dbh->prepare("SELECT * from gene where chromosome='X' and
	feature_type='GENE' and group_label='reference'");

#my $sqr1 = $dbh->prepare("SELECT * from repeat where chromosome='X' and feature_name like '%(%)n%'");
$sqr->execute();
open(WWW,">F:/Genius/STR/X_STR/Xgene.txt");
while(my $ref = $sqr->fetchrow_hashref())
{
    print WWW join("\t",$ref->{'chr_start'},$ref->{'chr_stop'},$ref->{'contig'},
    	$ref->{'contig_start'},$ref->{'contig_stop'},$ref->{'feature_name'},
        $ref->{'feature_id'},$ref->{'feature_type'}),"\n";
}
$dbh->disconnect();