#!/usr/bin/perl -w
use strict;

my $dir = "/home/orion/Documents/MT/seq/";
my $pmid = "15317881";
open(PAIR, $dir.$pmid."_pair") or die("Cannot open pair_match file");
open(MUT, $dir.$pmid) or die("Cannot open mutation file");
my %line_arr;
while(<MUT>) {
	chomp;
	if($_ ne "") {
		my $line = $_;
		#$line =~ s/[^_0-9a-zA-Z\.\t\|]//g;
		$line =~ /\b([A-Z]{1,3}\d+)\b/;
		$line_arr{$1} = $line;
	}
}
open(OUT, ">".$dir."tmp") or die("cannot write to outfile");
while(<PAIR>) {
	chomp;
	$_ =~ s/[^0-9a-zA-Z\t]//g;
	my @onepair = split("\t",$_);
	if(!exists($line_arr{$onepair[0]}) || !exists($line_arr{$onepair[1]})) {
		die("Error, Acc not found");
	}
	my @ele_arr1 = split("\t",$line_arr{$onepair[0]});
	my @ele_arr2 = split("\t",$line_arr{$onepair[1]});
	shift(@ele_arr1);
	shift(@ele_arr1);
	my $pos1 = shift(@ele_arr1);
	shift(@ele_arr2);
	shift(@ele_arr2);
	my $pos2 = shift(@ele_arr2);
	if($pos1 <= $pos2) {
		print OUT $onepair[0]."+".$onepair[1]."\t";
		print OUT join("\t",@ele_arr1,@ele_arr2);
	} else {
		print OUT $onepair[1]."+".$onepair[0]."\t";
		print OUT join("\t",@ele_arr2,@ele_arr1);	
	}
	print OUT "\n";
}