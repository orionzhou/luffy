#!/usr/bin/perl -w
require "seq_ret.pl";
use strict;
use Switch;
open(STS,"F:/Genius/STR/X_STR/sts_DXS_final.txt");
my $j=1;
my @sts_info;
my @seq_arr;
my $seq = "";
while(<STS>)
{
  chomp;
  @sts_info = split("\t",$_);
  @seq_arr = seq_ret($sts_info[1],$sts_info[2],23,0);
  $seq = $seq_arr[1];
  my $length = length($seq);
  my $cha = "";
  open(WWW,">F:/Genius/STR/DSP/".sprintf('%04.0f',$j).".txt");
  for(my $i=1;$i<=$length;$i++)
  {
    $cha = substr($seq,$i-1,1);
    switch ($cha)
    {
      case "A" {print WWW join("\t",1,0,0,0,"\n");}
      case "T" {print WWW join("\t",0,1,0,0,"\n");}
      case "C" {print WWW join("\t",0,0,1,0,"\n");}
      case "G" {print WWW join("\t",0,0,0,1,"\n");}
      else {print "Invalid Code!";last;}
    }
  }
  $j++;
}