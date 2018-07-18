#!/usr/bin/perl -w
use Switch;
use strict;
use Bio::Seq;
use Bio::SeqIO;
my $in = Bio::SeqIO->new(-file => "F:/Genius/STR/DSP/C.elegans.fasta", -format=>'fasta');
open(WWW,">F:/Genius/STR/DSP/seq_digital.txt");
while(my $seq = $in->next_seq())
{
  print "Sequence ",$seq->id,"(",$seq->length,")";
  my $buffer = $seq->subseq(7021,15020);
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
}