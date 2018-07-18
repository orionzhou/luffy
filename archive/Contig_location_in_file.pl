#!/usr/bin/perl -w
use strict;
my @chr_array = (1..22,"X","Y");
my $search_dir = "F:/Genius/STR/perl+php+fasta/human_genome_ref/";
foreach my $chr (@chr_array)
{
print "chromosome $chr :\n";
my $search_file = "ref_chr" . $chr . ".fa";
my $search_info = "CHR_" . $chr . ".txt";
my $search_file_path = $search_dir . $search_file;
my $search_info_path = $search_dir . $search_info;
open( S_FILE, $search_file_path ) or die("Could not open Search_File.");
open( S_INFO, $search_info_path ) or die("Could not open the AGP file");
open( W_INFO, ">".$search_dir."CHR".$chr.".txt" ) || die("Coundn't write file");
my $contig_count = 0;
while(<S_INFO>)
{
  chomp;
  if($_ =~ /^[0-9]+\t[0-9]+\t/)
  {
    my @contig_info = split( "\t" , $_ );
    my $part_print = ++$contig_count;
    my $contig_name = $contig_info[2];
    print join( "\t", "$part_print",
    	$contig_info[2], "$contig_info[0]-$contig_info[1]" ), "\n";
    seek(S_FILE, 0, 0);
    while(<S_FILE>)
    {
    	if($_ =~ /$contig_name/)
        {
          my $start_abs = tell(S_FILE);
          my $first_line = "";
          read(S_FILE, $first_line, 70);
          print "\tfrom $start_abs:$first_line\n";
          print W_INFO join("\t",$contig_info[0],$contig_info[1],$start_abs,$contig_info[2]),"\n";
          last;
        }
    }

  }
}
}