#!/usr/bin/perl -w
use Bio::Seq;
use Bio::SeqIO;
use warnings;
use strict;
my $search_dir = "F:/Genius/STR/perl+php+fasta/human_genome_ref/";
sub seq_ret()
{
	my ($pos_start, $pos_end, $chr_number, $flanking_len) = @_;
        if(!$flanking_len)
        {
        	$flanking_len = 0;
        }
        if( $pos_start !~ /^[0-9]{1,10}$/ || $pos_end !~ /^[0-9]{1,10}$/
        	 || $pos_start>=$pos_end )
        {
        	die("Invalid input parameter!");
        }
        elsif( abs($pos_start-$pos_end)<5 || abs($pos_start-$pos_end)>2000 )
        {
        	die("5-2000 bp allowed!");
        }
        else
        {
                my @chr_array = (1..22,'X','Y');
         	my $chr = $chr_array[$chr_number-1];
    		my $search_file = "ref_chr" . $chr . ".fa";
    		my $search_info = "CHR" . $chr . ".txt";
    		my $search_file_path = $search_dir . $search_file;
    		my $search_info_path = $search_dir . $search_info;
    		open( S_FILE, $search_file_path ) or die("Could not open Search_File.");
    		open( S_INFO, $search_info_path ) or die("Could not open the Info file");
                my $tag=-1;
                my $contig_start_rel, my $contig_end_rel, my $contig_start_abs, my $contig_name;
                while(<S_INFO>)
                {
                	chomp;
                        my @contig_info = split("\t",$_);
                        if($contig_info[0]<=$pos_start && $contig_info[1]>=$pos_end)
                        {
                        	$contig_start_rel = $contig_info[0];
                                $contig_end_rel = $contig_info[1];
                                $contig_start_abs = $contig_info[2];
                                $contig_name = $contig_info[3];
                                $tag = 1;
                        }
                }
                if($tag == -1)
                {
                	die("Request Sequence is not available through NCBI");
                }
                else
                {
                	my $up_len = $flanking_len<=($pos_start-$contig_start_rel)? $flanking_len : ($pos_start-$contig_start_rel);
                        my $down_len = $flanking_len<=($contig_end_rel-$pos_end)? $flanking_len : ($contig_end_rel-$pos_end);
                        my $seq_len = $pos_end - $pos_start + 1;
                        my $gap = $pos_start - $up_len - $contig_start_rel;
                        my $segment_start = $contig_start_abs + (int(($gap)/70))*71 + ($gap)%70;
                        my $up_len_true = (int($up_len/70))*71 + $up_len%70;
                        my $down_len_true = (int($down_len/70))*71 + $down_len%70;
                        my $seq_len_true = (int($seq_len/70))*71 + $seq_len%70;
                        #print join("\t",$up_len,$seq_len,$down_len),"\n";
                        #print join("\t",$up_len_true,$seq_len_true,$down_len_true),"\n";
                        seek(S_FILE, $segment_start, 0);
                        my $up = "";
                        my $down = "";
                        my $seq = "";
                        my $all = "";
                        read(S_FILE, $all, $up_len_true+$seq_len_true+$down_len_true+30);
                        $all =~ s/\W//g;
                        $up = substr($all,0,$up_len);
                        $seq = substr($all,$up_len,$seq_len);
                        $down = substr($all,$up_len+$seq_len,$down_len);
                        my @seqs = ($up,$seq,$down);
                        return @seqs;
                }
    	}
}
