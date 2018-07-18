#!/usr/bin/perl -w
#given a sequence dir, generate an alignment result file
require "./al2seq3.pl";
use strict;
use Bio::Seq;
use Bio::SeqIO;

my $seq_dir = "./";
my $option_included = shift;  # 1:entire genome, 2:HVRI&II, 3:Coding Region, 4:HVRI, 5:HVRII, 6:ATP6&8, 7:12s rRNA
my $option_excluded = shift;
my $seq_path = shift;
if(!$seq_path) {
	$seq_path = "IN";
}
if(!$option_included) {
	$option_included = 1;
}
if(!$option_excluded) {
	$option_excluded = 1;
}
my %range_included = (
		1 => "1-16569", 
		2 => "40-370,16024-16365",
		3 => "576-16023",
		4 => "16024-16365",
		5 => "40-370",
		6 => "8350-9250",
		7 => "640-1600"
		);
my %range_excluded = (
		1 =>"^303-315"
		);
my $rCRS_path = "rCRS.gb";

my $rCRS_file = Bio::SeqIO->new(-file => $rCRS_path, -format=>'genbank');
my $rCRS = $rCRS_file->next_seq();
my $seq_file = Bio::SeqIO->new( -file=>$seq_dir.$seq_path, -format=>'fasta' )
	or die("Cannot open sequence file");
open(STDOUT, ">OUT") or die("cannot write to outfile");
while(my $seq = $seq_file->next_seq()) {
	if($seq->length > 16000) {
		if($option_included == 5) {
			$seq = $seq->trunc(1,400);
		} elsif ($option_included==4 || $option_included==6 || $option_included==7) {
			my @pos = split("-",$range_included{$option_included});
			$seq = $seq->trunc(($pos[0]-50>0)?$pos[0]-50:1,($pos[1]+50<$seq->length)?$pos[1]+50:$seq->length);
		}
	}
	&al2seq($seq, $rCRS, split(",",$range_included{$option_included}), split(",",$range_excluded{$option_excluded}));
}
