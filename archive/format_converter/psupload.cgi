#!/usr/bin/perl -w
require "converter.pl";
use strict;
use CGI;
use Switch;
## Start setting up options here:
## Your path to where you want your files uploaded.
## Note: NO trailing slash
my $basedir = "/var/www/html/webroot/personal/genius/app/tmp";
## Do you wish to allow all file types? yes/no (no capital letters)
my $allowall = "no";
## If the above = "no"; then which is the only extention to allow?
## Remember to have the LAST 4 characters i.e. .ext
my $theext = ".txt";
## The page you wish it to forward to when done:
## I.E. http://www.mydomainname.com/thankyou.html
my $donepage = "http://www.xjtu.edu.cn/";
################################################
## DO NOT EDIT OR COPY BELOW THIS LINE ##
################################################
my $obj = new CGI;
print $obj->header(-type=>"text/html; charset=gb2312") ;
#print $obj->start_html( -title=>"Retrieving_Sequence", -BGCOLOR=>"#DFFCFD");
#print "Content-type: text/html\n";
my $onnum = 1;
my @fname_arr;
while ($onnum < 6)
{
	my $file = $obj->param("FILE$onnum");
	if ($file)
        {
		my $fileName = $file;
		$fileName =~ s/^.*(\\|\/)//;
                my $newmain = $fileName;
                my $filenotgood;
		if ($allowall ne "yes")
                {
		  if (lc(substr($newmain,length($newmain) - 4,4)) eq $theext)
                  {
		    $filenotgood = "yes";
		  }
		}
		if ($filenotgood eq "yes")
        	{
                  $fileName = join("_",$onnum,$fileName);
		  open (OUTFILE, ">$basedir/$fileName") || DIE("couldn't create file");
                  binmode(OUTFILE);
		  while (my $bytesread = read($file, my $buffer, 1024))
                  {
		    print OUTFILE $buffer;
		  }
		  close (OUTFILE);
                  push (@fname_arr,$fileName);
                  print "$fileName format ok, upload complete.<br>";
		}
                else
                {
                  print "Invalid File Format! Discarded!<br>";
                }
	}
        elsif ($onnum == 1)
        {
        	print "No file selected!<BR>";
                exit;
        }
	$onnum++;
}

my $filetype = $obj->param("fileType");
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
  #print "$fname\n";
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
      my $outname = substr($fname,0,length($fname)-4).".arp";
      push (@output, $outname);
    }
    case 2
    {
      ToGDA($basedir,$fname,$locus_count,$loci,$pcount,@genotype);
      my $outname = substr($fname,0,length($fname)-4).".nex";
      push (@output, $outname);
    }
    else { print "Invalid Outfile Format!"; }
  }
}

print $obj->h2({-align=>"left"}, "right click and save as...");
print "<p align='left'>";
foreach my $outfile (@output)
{
  print "<a href='./tmp/".$outfile."'>$outfile</a><br>";
}
print "</p>";