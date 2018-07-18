#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Net::SSH::Perl;
use Common;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;

my @KEYFILE = ("$ENV{'HOME'}/.ssh/id_rsa");
my $ssh = Net::SSH::Perl->new("lab.msi.umn.edu",
  debug => 1,
  identify_files => \@KEYFILE);
$ssh->login('zhoup');
my ($stdout, $stderr, $exit) = $ssh->cmd("pwd");
print $stdout;

