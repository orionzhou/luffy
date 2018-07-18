#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Switch;

my $working_dir = "./";
my $rCRS_path = Bio::SeqIO->new(-file => '/home/orion/Scripts/test/align/rCRS.gb', -format => 'genbank');
my $in_path = $working_dir."in";
my $out_path = $working_dir."out";
my $rCRS = $rCRS_path->next_seq();

print "=============================================\n"
	."1:\tgi|13272682|gb|AF346972.1|\t588\t1648\t750G\t1041G\t1438G\n"
	."2:\t(Sample ID\tExpected HG\tEstimated HG\tnp 16024-16569	np 001-437	np 438-576)\n"
	."\tCHN.ASN.0046\tN9a\tM33b\t16173C 16223\t73 150\t\n"
	."3:\t(Sample No.\tHG\tnp 16024-16569\tnp 001-437\tnp 438-576)\n"
	."\t8\tM7b1\t16129A	16192T\t73G 150T -249\t489C\n"
	."=============================================\n"
	."select input file format(default 1):";
my $input_format = <STDIN>;
chomp($input_format);
$input_format = ($input_format eq "") ? 1:$input_format; 
if($input_format !~ m/^[1-3]$/) {
	print "only options 1-3 are allowed";
	exit;
}
print "select output file format(default 2):";
my $output_format = <STDIN>;
chomp($output_format);
$output_format = ($output_format eq "") ? 2:$output_format;
if($output_format !~ m/^[1-3]$/) {
	print "only options 1-3 are allowed\n";
	exit;
}

my %range_include = (
	1 => "1-16569", 
	2 => "16024-16569,1-437,438-576",
	3 => "577-16023",
	4 => "16024-16569",
	5 => "1-437",
	6 => "438-576",
	7 => "8350-9250",
	8 => "640-1600"
	);
print "=============================================\n"
	."\t1(mtGenome):\t1-16569\n"
	."\t2(control region):\t16024-16569,1-437,438-576\n"
	."\t3(coding region):\t577-16023\n"
	."\t4(HVRI):\t16024-16569\n"
	."\t5(HVRII):\t1-437\n"
	."\t6(HVRIII):\t438-576\n"
	."\t7(ATP 6&8):\t8350-9250\n"
	."\t8(12s rRNA):\t640-1600\n"
	."=============================================\n"
	."select the range you want to include in output(default 1):";
my $range_include_option = <STDIN>;
chomp($range_include_option);
$range_include_option = ($range_include_option eq "") ? 1:$range_include_option;
if($range_include_option !~ m/^[1-8]$/) {
	print "only options 1-8 are allowed\n";
	exit;
}

my %range_exclude = (
	1 =>"303-315",
	2 =>""
	);
print "=============================================\n"
	."\t1(poly-C region):\t303-315\n"
	."\t2(none):\t\n"
	."=============================================\n"
	."select a range you want to exclude in output(default 1):";
my $range_exclude_option = <STDIN>;
chomp($range_exclude_option);
$range_exclude_option = ($range_exclude_option eq "") ? 1:$range_exclude_option;
if($range_exclude_option !~ /^[1-2]$/) {
	print "only options 1-2 are allowed\n";
	exit;
}

print "=============================================\n"
	."Want to exclude the ambiguous nucleotide substitutions, e.g. \"522N\"?\n"
	."=============================================\n"
	."y/n(default y):";
my $ignore_option = <STDIN>;
chomp ($ignore_option);
$ignore_option = ($ignore_option eq "") ? "y":$ignore_option;
if($ignore_option !~ /^[yn]$/) {
	print "only \"y\" or \"n\" are allowed\n";
	exit;
}
my @ignore_chrs;
if($ignore_option eq "y") {
	push (@ignore_chrs, "N");
	#push "p" when $input_format=2;
}

open( IN, "<".$in_path ) or die("Cannot open input file");
open( OUT, ">".$out_path ) or die("Cannot create output file");

switch ($input_format) {
	case 1 {
		&convert_1();
	} case 2 {
		&convert_2();
	} case 3 {
		&convert_3();
	} else {
		print "invalid input_format option\n";
		exit;
	}
}

sub convert_1() {
	&print_firstline();
	while(<IN>) {
		chomp;
		if($_) {
			my @ele_arr = split("\t",$_);
			my $id = shift(@ele_arr);
			shift(@ele_arr);
			shift(@ele_arr);
			&print_head($id);
			foreach my $ele (@ele_arr) {
				if($ele =~ /^([.0-9]+)([A-Zd])$/) {
					&check_write_ele($1,$2);
				} else {
					print "error string value on $id\n";
				}
			}
			print OUT "\n";
		}
	}
}

sub convert_2() {
	&print_firstline();
	while(<IN>) {
		chomp;
		if($_ && $_!~/^Sample/) {
			my @ele_arr = split("\t",$_);
			my $id = shift(@ele_arr);
			my $group = shift(@ele_arr);
			shift(@ele_arr);
			my %mut_arr;
			for(my $i=1;$i<=2;$i++) {
				my $tmp = shift(@ele_arr);
				if($tmp) {
					my @ele_arr2 = split(" ",$tmp);
					foreach my $ele2 (@ele_arr2) {
						if($ele2 =~ /^([.0-9]+)([A-Zd])$/) {
							$mut_arr{$1} = $2;
						} elsif($ele2=~/^(\d+)$/) {
							$mut_arr{$1} = Bio::Seq->new(-seq=>$rCRS->subseq($1,$1));
						} elsif($ele2=~/p/) {
							#do nothing
						} else {
							print "error string value on $id\n";
							exit;
						}
					}
				}
			}
			&print_head($id);
			foreach (sort {$a <=> $b } keys(%mut_arr)) {
				&check_write_ele($_,$mut_arr{$_});
			}
			print OUT "\n";
		}
	}
}

sub convert_3() {
	&print_firstline();
	while(<IN>) {
		chomp;
		if($_) {
			my @ele_arr = split("\t",$_);
			my $number = shift(@ele_arr);
			my $group = shift(@ele_arr);
			for(my $i=1;$i<=$number;$i++) {
				my %mut_arr;
				foreach my $ele (@ele_arr) {
					if($ele =~ /^([.0-9]+)([A-Z])$/) {
						$mut_arr{$1} = $2;
					} elsif($ele=~/^\-(\d+)$/) {
						$mut_arr{$1} = "d";
					} else {
						print "error string value on $number $group\n";
						exit;
					}
				}
				my $id = ($i==1)?$group:$group.".".$i;
				&print_head($id);
				foreach (sort {$a <=> $b } keys(%mut_arr)) {
					&check_write_ele($_,$mut_arr{$_});
				}
				print OUT "\n";
			}
		}
	}
}

sub check_write_ele() {
	my ($pos,$value) = @_;
	#check position
	my @range_included = split(",",$range_include{$range_include_option});
	my @range_excluded = split(",",$range_exclude{$range_exclude_option});
	my @include_begin;
    my @include_end;
    if(@range_included) {
		foreach my $ele (@range_included) {
			my @pair = split("-",$ele);
			push (@include_begin,$pair[0]);
			push (@include_end,$pair[1]);
		}
	} else {
		push (@include_begin,1);
		push (@include_end,$rCRS->length);
	}
	my @exclude_begin;
	my @exclude_end;
	if(@range_excluded) {
		foreach my $ele (@range_excluded) {
			my @pair = split("-",$ele);
			push (@exclude_begin,$pair[0]);
			push (@exclude_end,$pair[1]);
		}
	} else {
		push (@exclude_begin,0);
		push (@exclude_end,0);
	}
	my $flag_pos = 0;
	my $i = 0;
	foreach (@include_begin) {
		if($pos>=$_ && $pos<=$include_end[$i]) {
			$flag_pos = 1;
		}
		$i++;
	}
	my $j = 0;
	foreach (@exclude_begin) {
		if(int($pos)>=$_ && int($pos)<=$exclude_end[$j]) {
			$flag_pos = 0;
		}
		$j++;
	}
	
	my $flag_value = 1;
	foreach my $ignore_chr (@ignore_chrs) {
		if($value =~ /$ignore_chr/) {
			$flag_value = 0;
		}
	}
	
	if($flag_pos==1 && $flag_value==1) {
		switch ($output_format) {
			case 1 {
				&printw_1($pos,$value);
			} case 2 {
				&printw_2($pos,$value);
			} case 3 {
				&printw_3($pos,$value);
			} else {
				print "invalid output_format option\n";
				exit;
			}
		}
	}
}

sub print_firstline() {
	switch ($output_format) {
		case 1 {
		} case 2 {
			print OUT join("\t","Sample","Description","np 16024-576"),"\n";
		} case 3 {
			print OUT join("\t","Sample No.","Sample Id","np 16024-576"),"\n";
		} else {
			print "invalid output_format option\n";
			exit;
		}
	}
}

sub print_head() {
	my ($id) = @_;
	switch ($output_format) {
		case 1 {
			my @tmp_arr = split("-",$range_include{$range_include_option});
			if($range_include_option==2) {
				@tmp_arr=(16024,576);
			}
			if(@tmp_arr!=2) {
				print "Error in determining output range\n";
				exit;
			}
			print OUT join("\t",$id,@tmp_arr),"\t";
		} case 2 {
			print OUT join("\t",$id,""),"\t";
		} case 3 {
			print OUT join("\t",1,$id),"\t";
		} else {
			print "invalid output_format option\n";
			exit;
		}
	}
}

sub printw_1() {
	my ($pos,$value) = @_;
	print OUT $pos.$value."\t";
}

sub printw_2() {
	my ($pos,$value) = @_;
	print OUT $pos.$value." ";
}

sub printw_3() {
	my ($pos,$value) = @_;
	if($value eq "d") {
		print OUT "-".$pos." ";
	} else {
		print OUT $pos.$value." ";
	}
}