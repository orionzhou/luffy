#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  netinfo.pl - extract level/type/parent info from a NET file

=head1 SYNOPSIS
  
  netinfo.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Net)
    -o (--out)    output file (Tbl)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Location;

my ($fi, $fo) = ('') x 2;
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

print $fho join("\t", qw/id par lev type/)."\n";

my $hi = {};
my $ph = {0=>""};
while( <$fhi> ) {
  chomp;
  next if /(^\#)|(^\s*$)/;
  if(/(^ +)(fill.*)$/) {
    my $lev = (length($1) + 1) / 2;
    my @ps = split(" ", $2);
    my $n = (@ps - 7) / 2;
    my %h = map {$ps[7+$_*2] => $ps[7+$_*2+1]} (0..$n-1);
    my $id = $h{"id"};
    my $type = $h{"type"};
    exists $ph->{$lev-1} || die "no level $lev-1\n";
    my $pid = $ph->{$lev-1};

    $hi->{$id} ||= 0;
    $hi->{$id} ++;
    $id .= ".".$hi->{$id} if $hi->{$id} > 1;

    print $fho join("\t", $id, $pid, $lev, $type)."\n";

    $ph->{$lev} = $id;
  }
}
close $fhi;
close $fho;


__END__
