#!/usr/bin/perl
use DBI;
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy", {'RaiseError' => 1});

 my $query = "select * from SSR where chr='X' and feature_name REGEXP '^\\\\([A-Z]{3,4}\\\\)n\$'";
    my $sqr = $dbh->prepare($query);
    $sqr->execute();
    while(my $ref = $sqr->fetchrow_hashref())
    {
      my $ssr_start = $ref->{'chr_start'};
      my $ssr_stop = $ref->{'chr_stop'};
      my $ssr_id = $ref->{'id'};
      my $ssr_link = "http://202.117.161.53/app/contents/show_ssr.php?id=".$ssr_id;
      my $repeat_unit = length($ref->{'feature_name'}) - 3;
      my $repeat_num = int(($ssr_stop - $ssr_start) / $repeat_unit);
      print "\t";
      print join("\t",$ssr_id,$ref->{'feature_name'},$ssr_start,$ssr_stop,$repeat_unit,$repeat_num,$ssr_link);
      print "\n";
    }
