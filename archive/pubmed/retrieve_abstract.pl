#!/usr/bin/perl -w
use LWP::Simple;
use strict;

my $eutils = "http://www.ncbi.nlm.nih.gov/entrez/eutils/";
my $paper_dir = "E:/Genius/Papers";
my $script_dir = "E:/Scripts/pubmed";

&give_me_id($paper_dir, 0);
&batch_fetch();

sub batch_fetch
{
  opendir( IDDIR, $script_dir."/ids" ) or die("Could not open ids_dir\n");
  my $entity;
  while( $entity = readdir(IDDIR) )
  {
    if( $entity ne "." && $entity ne ".." )
    {
      open( IDFILE, $script_dir."/ids/".$entity ) or die("Can't open id_file\n");
      open( INFOFILE, ">".$script_dir."/info/".$entity ) or die("Can't create info_file\n");

      print $entity,"\n";
      while( !eof(IDFILE) )
      {
        my @id_arr = ();
        for(my $i=0; $i<5; $i++)
        {
          my $buffer = <IDFILE>;
          if( defined($buffer) )
          {
            chomp($buffer);
            push @id_arr, $buffer;
          }
        }
        &fetch_by_id( join(",",@id_arr) );
      }
    }
  }
}

sub fetch_by_id
{
  my ($id_string) = @_;
  my $efetch = "efetch.fcgi?db=pubmed"
  	."&id=".$id_string."&retmode=text&rettype=medline";
  my $result = get($eutils . $efetch);
  print INFOFILE $result;
  print "\t5 done\n";
}

sub give_me_id()
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
      if( $entity =~ /^([0-9]{6,9}).pdf$/ )
      {
        my $id_file = $script_dir."/ids/".join(".",split("/",substr($dir,17)));
        print "\t" x $level,"|___",$entity,"\n";
        if( -e $id_file.".txt" )
        {
          open( ID, ">>".$id_file.".txt" )
        	or die("Could not open id_list file!\n");
        }
        else
        {
          open( ID, ">".$id_file.".txt" )
        	or die("Could not open id_list file!\n");
        }
        print ID $1,"\n";
        close( ID );
      }
      elsif( $entity eq "id.txt" || $entity eq "info.txt" )
      {
        unlink( $dir."/".$entity );
      }
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
      print "\t" x $level,"|_",$dir_to_access,"\n";
      &printdir( $dir."/".$dir_to_access, $level+1 );
    }
  }
}