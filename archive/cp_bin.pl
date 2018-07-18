#!/usr/bin/perl -w
use strict;
my $dir_m = $ENV{"m"};
my $dir_o = $ENV{"h"}."/bin";
print "copying files from $dir_m to:\n";
print "\t$dir_o\n";

my @fs;
push @fs, get_binaries("$dir_m/mt3");
push @fs, get_binaries("$dir_m/mt4");

my $f_str = join(" ", @fs);
system("cp $f_str $dir_o");

sub get_binaries {
  my ($dir) = @_;
  my @fs;
  opendir(my $dh, $dir) || die "cannot open $dir\n";
  for (readdir($dh)) {
    my $f = "$dir/$_";
    push @fs, $f if -f $f && -B $f;
  }
  return @fs;
}

