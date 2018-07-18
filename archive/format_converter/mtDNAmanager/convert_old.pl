#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;

my $working_dir = "./";
my $ref = Bio::SeqIO->new(-file => '/home/orion/Scripts/test/align/rCRS.gb', -format => 'genbank');
my $in_file = $working_dir."in";
my $out_file = $working_dir."out";

my $rCRS = $ref->next_seq();
#my $HVR1 = $rCRS->subseq(16024, 16384);

open( IN, "<".$in_file ) or die("Cannot open input file");
open( OUT, ">".$out_file ) or die("Cannot create output file");

my $line = <IN>;
print OUT $line;
while( $line = <IN> ) {
	chomp($line);
	my @ele_arr = split( "\t", $line );
	my $key=0;
	foreach my $value (@ele_arr) {
		if($key==0) {
			print OUT $value."\t";
		} elsif($key==1 || $key==2) {
			print OUT $value."\t";
		} elsif($key==3) {
			my @mut_arr_ori = split(" ",$value);
			foreach my $mut (@mut_arr_ori) {
				print OUT &convert($mut)."\t";
			}
		} elsif($key==4) {
			print OUT "\n\t\t\t";
			my @mut_arr_ori = split(" ",$value);
			foreach my $mut (@mut_arr_ori) {
				print OUT &convert($mut)."\t";
			}
		} elsif($key==5) {
			print OUT "\n\t\t\t";
			my @mut_arr_ori = split(" ",$value);
			foreach my $mut (@mut_arr_ori) {
				print OUT &convert($mut)."\t";
			}
		}
		$key ++;
	}
	print OUT "\n";
}
sub convert() {
	my ($string) = @_;
	if($string =~ /^(\d+)([dA-Z]{0,1})([.]{0,1})([p1-9]{0,1})([A-Z]{0,1})$/) {
		if($1.$2.$3.$4.$5 ne $string) {
			return "Error";
		} else {
			if($3 eq "") {  #not insertion
				if($2 eq "d") { #deletion
					return $1."[".$rCRS->subseq($1,$1)."/-]";
				} elsif ($2 eq "") { #transition
					return $1."[".$rCRS->subseq($1,$1)."/".&seq_comp($rCRS->subseq($1,$1))."]";
				} else { #transverstion
					return $1."[".$rCRS->subseq($1,$1)."/".$2."]";
				}
			} else { #insertion
				if($4 eq "p") {
					return $1."^".($1+1)."[-/poly".$5."]";
				} else {
					return $1."^".($1+1)."[-/".&char_gen($5,$4)."]"
				}
			}
		}
	} else {
		return "Error";
	}
}

sub seq_comp() {
	my ($nt) = @_;
	my $seq = Bio::Seq->new( -display_id => 'Original_NT', -seq => $nt);
	return $seq->revcom()->seq();
}

sub char_gen() {
	my ($str, $times) = @_;
	my $rs = "";
	for(my $i=0; $i<$times; $i++) {
		$rs .= $str;
	}
	return $rs;
}
