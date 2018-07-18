#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;

my $working_dir = "./";
my $rCRS_path = Bio::SeqIO->new(-file => '/home/orion/Scripts/align/rCRS.gb', -format => 'genbank');
my $in_path = $working_dir."in";
my $out_path = $working_dir."out";
my $rCRS = $rCRS_path->next_seq();

open( IN, "<".$in_path ) or die("Cannot open input file");
open( OUT, ">".$out_path ) or die("Cannot create output file");

my $line_count = 0;
my @id_arr;
while(<IN>) {
	chomp;
	last if !$_;
	$line_count ++;
	if($line_count == 1) {
		next;
	} elsif($line_count == 2) {
		@id_arr = split("\t",$_);
		splice(@id_arr,0,3);
		next;
	}
	my @ele_arr = split("\t",$_);
	my $pos = $ele_arr[1];
	my $value = $ele_arr[2];
	if($value ne $rCRS->subseq($pos,$pos)) {
		print join("\t",$pos,$value,$rCRS->subseq($pos,$pos))."\n"; 
	}
}

for(my $sample_count=1; $sample_count<=@id_arr; $sample_count++) {
	print OUT $id_arr[$sample_count-1]."\t";
	seek(IN,0,0);
	my $line_count = 0;
	while(<IN>) {
		chomp;
		last if !$_;
		$line_count ++;
		if($line_count==1 || $line_count==2) {
			next;
		}
		my @ele_arr = split("\t",$_);
		my $pos = $ele_arr[1];
		my $rCRS_value = $ele_arr[2];
		my $snp = $ele_arr[$sample_count+2];
		if($rCRS_value ne $snp) {
			if($snp eq "DEL") {
				print OUT $pos."d ";
			} elsif(length($snp)>length($rCRS_value)) {
				if($rCRS_value ne substr($snp,0,1)) {
					print OUT $pos.substr($snp,0,1)." ";
					print "Mutation+Insertion found at $id_arr[$sample_count-1] - $pos\n";
				}
				for(my $i=1;$i<length($snp);$i++) {
					print OUT $pos.".".$i.substr($snp,$i,1)." ";
				}
			} else {
				print OUT $pos.$snp." ";
			}
		}
	}
	print OUT "\n";
}