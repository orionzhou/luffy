#!/usr/bin/perl -w
#use the alignment result file, generate an formatted sequence file (for arlequin use)
use strict;
use Bio::Seq;
use Bio::SeqIO;
my $in = Bio::SeqIO->new(-file => 'E:/Scripts/test/align/rCRS.gb', -format => 'genbank');
my $align_result_file = "E:/Scripts/test/align/align_result1.txt";

my $rCRS = $in->next_seq();
#my $HVR1 = $rCRS->subseq(16024, 16384);

open( ALIGN, $align_result_file ) or die("Cannot open align result file");

my $start = 16000; my $stop = 16400;
my @mut_arr_ori;
while( my $line = <ALIGN> ) {
	chomp($line);
	my @ele_arr = split( "\t", $line );
	my $key=0;
	foreach my $value (@ele_arr) {
		if($key == 0) {
			#print $value,"\t";
		} elsif($key == 1) {
			$start = ($start>$value) ? $start : $value;
		} elsif($key == 2) {
			$stop = ($stop<$value) ? $stop : $value;
		} else {
			if($value =~ /^(\d+)\^(\d+)\[\-\/([A-Z]+)\]$/) {
				push @mut_arr_ori, "$1^".$2."[$3]\n";
			}
		}
		$key ++;
	}
}
print "shortest segment: $start-$stop\n";

my %mut_arr = &clean(@mut_arr_ori);
seek(ALIGN, 0, 0);
while( my $line = <ALIGN> ) {
	chomp($line);
	if($line ne "") {
		&seq_recover(16001,16400,$line,%mut_arr);
		print "\n";
	}
}

sub clean() {
	my @mut_arr_ori = @_;
	my %mut_arr;
	foreach my $mut_ori (@mut_arr_ori) {
		$mut_ori =~ /^(.*)\^.*\[(.*)\]$/;
		#print $1."-".$2."(".length($2).")\n";
		$mut_arr{$1} = exists($mut_arr{$1}) ?
			(($mut_arr{$1}<length($2))?length($2):$mut_arr{$1}) : length($2);
	}
	foreach (sort {$a <=> $b } keys(%mut_arr)) {
	}
	return %mut_arr;
}

sub seq_recover() {
	my ($start, $stop, $line, %mut_arr_global) = @_;
	my @ele_arr = split("\t",$line);
	my %mut_arr;
	my $key = 0;

	print join("\t",@ele_arr[0..2]),"\t";
	foreach my $value (@ele_arr) {
		if($key >= 3) {
			if($value =~ /^(\d+)\^(\d+)\[\-\/([A-Z]+)\]$/) {
				while(my ($gg,$mm) = each(%mut_arr_global) ) {
					if($gg == $1) {
						my $tmp = exists($mut_arr{$gg}) ? $mut_arr{$gg} : $rCRS->subseq($gg,$gg);
						$mut_arr{$gg} = $tmp.$3.&char_gen("-",$mm-length($3));
						delete $mut_arr_global{$gg};
					}
				}
			} elsif($value =~ /^(\d+)\[(.)\/(.)\]$/) {
				if(exists($mut_arr{$1})) {
					$mut_arr{$1} = $3.substr($mut_arr{$1},1,length($mut_arr{$1}));
				} else {
					$mut_arr{$1} = $3;
				}
			} else {
				print "Unknown Mutation\n";
				exit;
			}
		}
		$key ++;
	}
	while(my ($gg,$mm) = each(%mut_arr_global) ) {
		my $tmp = exists($mut_arr{$gg}) ? $mut_arr{$gg} : $rCRS->subseq($gg,$gg);
		$mut_arr{$gg} = $tmp.&char_gen("-",$mm);
	}

	my $seqstr = "";
	my $start_loc = ($start<$ele_arr[1]) ? $ele_arr[1] : $start;
	my $stop_loc = ($stop>$ele_arr[2]) ? $ele_arr[2] : $stop;

	$seqstr .= &char_gen("?",$start_loc-$start);
	my $tmp_loc = $start_loc;
	foreach (sort {$a <=> $b } keys(%mut_arr)) {
		if($_ >= $start_loc && $_ <= $stop_loc) {
			#print $_."[".$mut_arr{$_}."]\t";
			$seqstr .= ($tmp_loc<$_) ? $rCRS->subseq($tmp_loc,$_-1) : "";
			$seqstr .= $mut_arr{$_};
			$tmp_loc = $_+1;
		} elsif( ($_>=$start && $_<$start_loc) || ($_>$stop_loc && $_<=$stop) ) {
			$seqstr .= &char_gen("?",length($mut_arr{$_})-1);
		}
	}
	$seqstr .= ($tmp_loc<$stop_loc) ? $rCRS->subseq($tmp_loc,$stop_loc) : "";
	$seqstr .= &char_gen("?",$stop-$stop_loc);

	#check with the origninal sequence
	#print "\n";
	#$seqstr =~ s/[\?\-]//g;
	#my $seqIO = Bio::SeqIO->new(-file => 'E:/Scripts/test/mito/mito_seq/'.$ele_arr[0].".gb",
	#	 -format => 'genbank');
	#my $seq = $seqIO->next_seq();
	#if($seq->seq() ne $seqstr) {
	#	if($seq->seq() !~ /$seqstr/) {
	#		print $seq->seq()."\n$seqstr\n";
	#	}
	#}
	print $seqstr;
}

sub char_gen() {
	my ($str, $times) = @_;
	my $rs = "";
	for(my $i=0; $i<$times; $i++) {
		$rs .= $str;
}
	return $rs;
}