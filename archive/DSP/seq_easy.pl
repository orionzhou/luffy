#!/usr/bin/perl -w
use strict;
use Switch;
open(WWW,">F:/Genius/STR/DSP/seq12_digital.txt");
open(IN,"F:/Genius/STR/DSP/seq12.txt");
my $buffer = <IN>;
chomp;
my $length = length($buffer);
my $cha = "";
for(my $i=1;$i<=$length;$i++)
{
    $cha = substr($buffer,$i-1,1);
    switch ($cha)
    {
      case "A" {print WWW join("\t",1,0,0,0,"\n");}
      case "T" {print WWW join("\t",0,1,0,0,"\n");}
      case "C" {print WWW join("\t",0,0,1,0,"\n");}
      case "G" {print WWW join("\t",0,0,0,1,"\n");}
      else {print "Invalid Code!";last;}
    }
}