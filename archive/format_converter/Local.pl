#!/usr/bin/perl -w
require "converter.pl";
use strict;
use Switch;
my $basedir = "D:/";
my @fname_arr = ( "9STR_achang.txt", "9STR_baoan.txt", "16STR_dulong.txt",
	"16STR_bai.txt", "9STR_chaoxian.txt", "9STR_dai.txt",
         "9STR_deang.txt", "9STR_dongxiang.txt", "9STR_ewenki.txt",
          "9STR_han.txt", "9STR_hazak.txt", "9STR_hui.txt",
           "9STR_jingpo.txt", "9STR_kerkezi.txt", "9STR_lisu.txt",
            "9STR_menggu.txt", "9STR_miao.txt", "9STR_naxi.txt",
             "9STR_nu.txt", "9STR_pumi.txt", "9STR_sala.txt",
              "9STR_tu.txt", "9STR_uzbik.txt", "9STR_wei.txt",
               "9STR_xibo.txt", "9STR_yao.txt", "9STR_yi.txt",
                "9STR_yugu.txt", "9STR_zang.txt", "9STR_zhuang.txt",
                 );
my $filetype = 1;

#my $filetype = $obj->param("fileType");
my @output;
foreach my $fname (@fname_arr)
{
  open(FI,"$basedir/$fname") || die("couldn't find upload file!");
  my $firstline = <FI>;
  if( $firstline !~ /^sample.*/i )
  {
    seek(FI,0,0);
  }
  my $pid = "";
  my $locus_count = 0;
  my @locus_name;
  print "$fname\n";
  while( <FI> )
  {
    chomp;
    my @ele_arr = split("\t",$_);
    if( scalar(@ele_arr)>0 and $ele_arr[1] !~ /amel/i )
    {
      if( $pid eq "" || $pid eq $ele_arr[0])
      {
        $pid = $ele_arr[0];
        $locus_name[$locus_count] = $ele_arr[1];
        $locus_count ++;
      }
      else
      {
        seek(FI,0,0);
        $firstline = <FI>;
        if( $firstline !~ /^sample.*/i )
        {
          seek(FI,0,0);
        }
        last;
      }
    }
  }
  #print "$locus_count loci: ",join(" ",@locus_name),"\n";
  $pid = "";
  my @genotype;
  my $pcount = 0;
  my $lcount = 0;
  my $g1 = "";
  my $g2 = "";
  while( <FI> )
  {
    chomp;
    my @ele_arr = split("\t",$_);
    switch(scalar(@ele_arr))
    {
      case 1 {$g1 = $g2 = "?";}
      case 2 {$g1 = $g2 = "?";}
      case 3
      {
        $g2 = "?";
        if($ele_arr[2] eq "" | $ele_arr[2] eq "0")
        {
          $g1 = "?";
        }
        else
        {
          $g1 = uc($ele_arr[2]);
        }
      }
      case 4
      {
        if($ele_arr[2] eq "" | $ele_arr[2] eq "0")
        {
          $g1 = "?";
        }
        else
        {
          $g1 = uc($ele_arr[2]);
        }
        if($ele_arr[3] eq "" | $ele_arr[3] eq "0")
        {
          $g2 = "?";
        }
        else
        {
          $g2 = uc($ele_arr[3]);
        }
      }
    }
    if( scalar(@ele_arr)>0 and $ele_arr[1] !~ /amel/i )
    {
      if( $pid eq "" || $pid ne $ele_arr[0])
      {
        $lcount = 0;
        $pid = $ele_arr[0];
        $pcount ++;
      }
      $lcount ++;
      $genotype[$pcount-1][$lcount-1] = join("*",$g1,$g2);
    }
  }
  #print "\tTotally $pcount samples.\n";
  my $loci = join("\t",@locus_name);
  switch ($filetype)
  {
    case 1
    {
      ToArn($basedir,$fname,$locus_count,$loci,$pcount,@genotype);
      push (@output, substr($fname,0,length($fname)-4).".arp");
    }
    case 2
    {
      ToGDA($basedir,$fname,$locus_count,$loci,$pcount,@genotype);
      push (@output, substr($fname,0,length($fname)-4).".nex");
    }
    else { print "Invalid Outfile Format!"; }
  }
}
print join("\n\t",@output);