#!/usr/bin/perl -w
open(RC,"F:/Genius/STR/X_STR/Xrepeat.txt");
open(Re,">F:/Genius/STR/X_STR/Xrepeat_ref.txt");
my $buffer = '';
while(<RC>)
{
	$buffer = $_;
	if($buffer =~ /reference/)
        {
        	print Re $buffer;
        }
}

