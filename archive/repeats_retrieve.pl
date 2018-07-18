#!/usr/bin/perl -w
require "seq_ret.pl";
open(RAW,"F:/Genius/STR/X_STR/Xrepeat_ref.txt");
open(WWW,">F:/Genius/STR/X_STR/Xrepeat_filtered.txt");
while(<RAW>)
{
	chomp;
        my @str_info = split("\t",$_);
        if( $str_info[9] =~ /\(([atcg]{3,4})\)n/i )
        {
        	my @seq_array = seq_ret($str_info[2],$str_info[3],23,0);
        	my $seq = $seq_array[1];
                my $seq_pattern = &seq2pattern($seq,$1);
                my $repeat_num = int(length($seq)/length($1));
		print WWW join("\t",@str_info[2,3],@str_info[5..7],
                	@str_info[9,10],$seq_pattern,$repeat_num,$seq,"\n");
        }
}
sub seq2pattern()
{
	my ($seq,$p1) = @_;
        my $len1 = length($p1);
        my $pattern = "";
        if( $seq =~ /(($p1){2,})/i )
        {
        	my $len2 = length($1);
                my $unit = $len2/$len1;
                $pattern = "(".$p1.")".$unit;
        }
        return $pattern;
}