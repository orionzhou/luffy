#!/usr/bin/perl -w
#given a sequence file, generate an alignment result file
require "./al2seq3.pl";
use strict;
use Bio::Seq;
use Bio::SeqIO;

my $seq_dir = "/home/orion/Documents/MT/seq/fasta/";
my $list_path = $seq_dir."list";
my $rCRS_path = "rCRS.gb";
my $pmid = shift;

open(LIST, $list_path) || die("cannot open list file");
my $rCRS_file = Bio::SeqIO->new(-file => $rCRS_path, -format=>'genbank');
my $rCRS = $rCRS_file->next_seq();

my $out_file = $seq_dir."../align_result_tmp";
open(STDOUT, ">".$out_file) || die("cannot write to outfile");
while(<LIST>) {
  chomp;
  my $line = $_;
  if($line =~ /$pmid/) {
    my @ele_arr = split("\t",$line);
    my $seq_file_path = $seq_dir;
    foreach my $ele (@ele_arr[2..(@ele_arr-1)]) {
      if($ele ne "") {
        if($ele !~ /\[.*\]/) {
          $seq_file_path .= "[$ele]";
        } else {
          $seq_file_path .= $ele;
        }
      }
    }
    $seq_file_path .= "[$ele_arr[0]-$ele_arr[1]].Fasta";
    my $seq_file = Bio::SeqIO->new( -file=>$seq_file_path, -format=>'fasta' ) || die("cannot open seq_file {$seq_file_path}.Fasata");
    while( my $seq = $seq_file->next_seq() ) {
      #my $seq_truncated = $seq->trunc(16001,$seq->length);
      my @range = ("16024-16365","40-302","316-370");
      &al2seq($seq, $rCRS, @range);
    }
  }
}
print "\n";
close(STDOUT);
#HVRI: 16001,$seq->length
#HVRII: 1,600
#ATP6/8: 8301,9300
