#!/usr/bin/perl -w
use strict;
open(STS,"F:/Genius/STR/X_STR/sts_DXS_final.txt");
open(REP,"F:/Genius/STR/X_STR/Xrepeat_final.txt");
open(STSREP,">F:/Genius/STR/X_STR/sts_DXS_rep.txt");
my @sts_info;
my @repeat_info;
my @ele_arr;
my $buffer = "";
my $tag = 0;
my $hit = 0;
while(<REP>)
{
  chomp;
  @ele_arr = split('\t',$_);
  $buffer = join("\t",@ele_arr[0..2]);
  push @repeat_info, $buffer;
}
my @repeat_arr;
while(<STS>)
{
  chomp;
  @sts_info = split("\t",$_);
  print STSREP join("\t",@sts_info[0..2,6,7]);
  $tag = 0;
  foreach my $line (@repeat_info)
  {
    @repeat_arr = split("\t",$line);
    if($sts_info[1]<=$repeat_arr[1] && $sts_info[2]>=$repeat_arr[2])
    {
      print STSREP "\t",join("\t",@repeat_arr);
      $tag = 1;
    }
  }
  if($tag == 1)
  {
    $hit ++;
  }
  print STSREP "\n";
}
print "\nTotally $hit hits.";