#!/usr/bin/perl -w
use strict;
my $k = "kkk";
sub ToArn
{
  my ($basedir,$fname,$locus_count,$loci,$pcount,@genotype) = @_;
  my @locus_name = split("\t",$loci);
  open ( WR, ">$basedir/".substr($fname,0,length($fname)-4).".arp" );
  print WR "[Profile]\n";
  print WR "\tTitle=\"User's Input Data\"\n";
  print WR "\t\tNbSamples=1\n";
  print WR "\t\tDataType=MICROSAT\n";
  print WR "\t\tGenotypicData=1\n";
  print WR "\t\tGameticPhase=0\n";
  print WR "\t\tMissingData=\"?\"\n";
  print WR "\t\tLocusSeparator=\"\t\"\n";
  print WR "[Data]\n";
  print WR "\t[[Samples]]\n";
  print WR "\t#There are $locus_count loci: ",join(",",@locus_name),"\n";
  my $test;
  my $epcount = $pcount;
  for(my $i=0; $i<$pcount; $i++)
  {
    $test = 0;
    while( $genotype[$i][$test] eq "?*?" )
    {
      $test ++;
      if($test == $locus_count)
      {
        $epcount --;
        $genotype[$i][0] = "none";
        last;
      }
    }
  }
  print WR "SampleName=\"Unknown\"\n";
  print WR "SampleSize=$epcount\n";
  print WR "SampleData={\n";
  my $ei = 0;
  for(my $i=0; $i<$pcount; $i++)
  {
    if($genotype[$i][0] ne "none")
    {
      print WR "S",$ei+1,"\t1\t",;
      my @temp;
      for(my $j=0; $j<$locus_count; $j++)
      {
        my @ele_2 = split('\*',$genotype[$i][$j]);
        print WR $ele_2[0],"\t";
        $ele_2[1] =~ s/\s//;
        $temp[$j] = $ele_2[1];
      }
      print WR "\n\t\t";
      for(my $j=0; $j<$locus_count; $j++)
      {
        print WR $temp[$j],"\t";
      }
      $ei ++;
      print WR "\n";
    }
  }
  print WR "}";
}

sub ToGDA
{
  my ($basedir,$fname,$locus_count,$loci,$pcount,@genotype) = @_;
  my @locus_name = split("\t",$loci);
  open ( WR, ">$basedir/".substr($fname,0,length($fname)-4).".nex" );
  print WR "#nexus\n\n";
  print WR "begin gdadata;\n";
  print WR "\tdimensions nloci=$locus_count npops=1;\n";
  print WR "\tformat interleaved nolabels separator=/ missing=?;\n";
  print WR "\tlocusallelelabels\n";
  for(my $k=0; $k<scalar(@locus_name); $k++)
  {
    print WR "\t\t",$k+1," $locus_name[$k]";
    if($k != scalar(@locus_name) - 1)
    {
      print WR ",";
    }
    print WR "\n";
  }
  print WR "\t;\n";
  print WR "\tmatrix\n";
  my $test;
  my $epcount = $pcount;
  for(my $i=0; $i<$pcount; $i++)
  {
    $test = 0;
    while( $genotype[$i][$test] eq "?*?" )
    {
      $test ++;
      if($test == $locus_count)
      {
        $epcount --;
        $genotype[$i][0] = "none";
        last;
      }
    }
  }
  print WR "\t\tPopulation:\n";
  my $ei = 0;
  for(my $i=0; $i<$pcount; $i++)
  {
    if($genotype[$i][0] ne "none")
    {
      print WR "\t\t\tS",$ei+1,"\t",;
      my @temp;
      for(my $j=0; $j<$locus_count; $j++)
      {
        my $temp = $genotype[$i][$j];
        $temp =~ s/\*/\//;
        print WR $temp,"\t";
      }
      $ei ++;
      print WR "\n";
    }
  }
  print WR "\t;\n";
  print WR "end;\n";
}