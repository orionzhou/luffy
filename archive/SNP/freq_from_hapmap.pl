#!/usr/bin/perl -w
use strict;
use LWP::Simple;
my @snp_arr = ("rs1805055");

foreach my $rs (@snp_arr)
{
  my $url = "http://www.hapmap.org/cgi-perl/snp_details?"
  	."name=$rs&source=hapmap_B35";
  my $html = get($url);
  $html =~ //;
  print $1;
}