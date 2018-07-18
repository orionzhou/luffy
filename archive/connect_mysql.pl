#!/usr/bin/perl -w
use DBI;
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy", {'RaiseError' => 1});

# query
#my $sqr = $dbh->prepare("SELECT * from sts where chromosome='X' and
#	feature_name like '%DX%S%' and group_label='reference'");

my $sqr = $dbh->prepare("SELECT * from repeat where chromosome='X' and
	group_label='reference'");
$sqr->execute();
open(WWW,">F:/Genius/STR/X_STR/Xrepeat_ori.txt");
while(my $ref = $sqr->fetchrow_hashref())
{
    print WWW join("\t",$ref->{'chr_start'},$ref->{'chr_stop'},$ref->{'contig'},
    	$ref->{'contig_start'},$ref->{'contig_stop'},$ref->{'feature_name'},
        $ref->{'feature_id'}),"\n";
}
$dbh->disconnect();