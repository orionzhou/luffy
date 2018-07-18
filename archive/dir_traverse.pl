#!/usr/bin/perl -w
use strict;

my $dir = "E:/Scripts";
&printdir( $dir, 0 );

sub printdir()
{
  my ($dir,$level) = @_;
  opendir( SEQ, $dir ) or die("Could not open seq_dir\n");
  my $entity;
  my $dircount = 0;
  my @dir_arr = ();
  while( $entity = readdir(SEQ) )
  {
    if( -d $dir."/".$entity )
    {
      if( $entity ne "." && $entity ne ".." )
      {
        push @dir_arr, $entity;
      }
      $dircount ++;
    }
    elsif( -f $dir."/".$entity )
    {
      print "\t" x $level,"File\t-\t",$entity,"\n";
    }
  }
  if( $dircount == 2 )
  {
    return;
  }
  else
  {
    foreach my $dir_to_access (@dir_arr)
    {
      print "\t" x $level,"Directory\t-\t",$dir_to_access,"\n";
      &printdir( $dir."/".$dir_to_access, $level+1 );
    }
  }
}