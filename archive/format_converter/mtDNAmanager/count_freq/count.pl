#!/usr/bin/perl -w
use strict;

my $working_dir = "./";
my $in_file = $working_dir."in";
my $out_file = $working_dir."out";

open( IN, "<".$in_file ) or die("Cannot open input file");
open( OUT, ">".$out_file ) or die("Cannot create output file");

my $prev="";
my $count=0;
my @id_arr;
while(<IN>) {
	chomp($_);
	if($_) {
		my @ele_arr = split("\t",$_);
		my $id = shift(@ele_arr);
		my $current = shift(@ele_arr);
		if($current ne $prev) {
			if($prev ne "") {
				print OUT join("\t",$prev,scalar(@id_arr),join(" ",@id_arr),"\n");
			}
			@id_arr = ();
		}
		push (@id_arr,$id);
		$prev = $current;
	}
}
print OUT join("\t",$prev,scalar(@id_arr),join(" ",@id_arr),"\n");