#!/usr/bin/perl -w
use strict;
my $file_prefix = "E:/Genius/Gene/Sequence/";

#my @name_arr = ("FMR1.txt","DMD.txt");
#my @name_arr = ("ABCA1.txt","APOA2.txt","APOA5.txt","APOC3.txt","APOE.txt","FASN.txt","LIPE.txt");
my @name_arr = ("DTNBP1.txt","NRG1.txt","DAOA.txt","PRODH.txt","COMT.txt","DISC1.txt");

foreach my $name (@name_arr)
{
  open( NA, $file_prefix.$name ) or die("Couldn't Open File.");
  print substr($name,0,length($name)-4);;
  my $temp = <NA>;
  while($temp !~ /^Enzyme\s+Freq/)
  {
    $temp = <NA>;
  }
  my @enzyme_arr;
  my @freq_arr;
  my @sites_arr;
  my $sites;
  while(<NA>)
  {
    chomp;
    if( $_ =~ /^(.*):(.*)$/ )
    {
      my $buffer = $2;
      if( $1 =~ /^([a-zA-Z0-9]+)\s+(\d+)/ )
      {
        #print "$1:$2 cuts\n";
        push(@enzyme_arr,$1);
        push(@freq_arr,$2);
      }
      while($buffer =~ /\s(\d+)\s/g)
      {
        #print "$1\t";
        $sites .= $1." ";
      }
    }
    else
    {
      #print "\n\n";
      push(@sites_arr,$sites);
      $sites = "";
    }
  }
  print "\t",join("\t",@enzyme_arr),"\n";
  print "\t",join("\t",@freq_arr),"\n";
  my @ele_count;
  foreach my $string (@sites_arr)
  {
    my @ele_arr = split(" ",$string);
    push(@ele_count,scalar(@ele_arr));
  }
  print "\t",join("\t",@ele_count),"\n";
  #print "\t",join("\t",@sites_arr),"\n";
}
print "\n";