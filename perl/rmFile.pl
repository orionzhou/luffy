#!/usr/bin/perl -w
use strict;
my $dir = "/project/youngn/zhoup/Scripts";
rm_hidden_files($dir);
sub rm_hidden_files {
    my ($dir) = @_;
    die "$dir not a directory" unless -d $dir;
    my @fps;
    opendir(JJ, $dir) or die "cannot open $dir\n";
    my @fns = readdir(JJ);
    for my $fn (@fns) {
        next if $fn eq "." || $fn eq "..";
        my $fp = "$dir/$fn";
        if( -f $fp ) {
            if($fn =~ /^\.\_/) {
                system("rm $fp");
                print "delete $fp\n";
            } else {
                push @fps, $fp;
            }
        } elsif( -d $fp) {
            if($fn eq ".svn") {
                system("rm -rf $fp");
                print "delete $fp\n";
            } else {
                push @fps, rm_hidden_files($fp);
            }
        }
    }
    return @fps;
}




