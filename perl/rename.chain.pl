#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  rename.chain.pl - rename fields with a map file for a Chain file

=head1 SYNOPSIS
  
  rename.chain.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -m (--map)    ID mapping file
    -p (--opt)    option (qry|tgt, default: tgt)

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
use List::Util qw/min max sum/;

my ($fi, $fo, $fm) = ('', '', '');
my $opt = "tgt";
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "map|m=s"  => \$fm,
  "opt|p=s"  => \$opt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fm;

$opt = lc($opt);
die "unknown opt: $opt\n" unless $opt =~ /^(qry)|(tgt)$/;

if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $h;
open(my $fhm, "<$fm") or die "cannot read $fm\n";
while(<$fhm>) {
  chomp;
  my @ps = split "\t";
  my ($oid, $nid) = @ps[0, -1];
  exists $h->{$oid} && die "$oid has >1 mappings\n";
  $h->{$oid} = $nid;
}
close $fhm;

while( <$fhi> ) {
  chomp;
  my $line = $_;
  if($line =~ /(^id)|(^\#)|(^\s*$)/) {
    print $fho $line."\n";
    next;
  }
  my @ps = split(" ", $line, -1);
  if(@ps < 13) {
    print $fho "$line\n";
    next;
  }
  
  if($opt eq "qry") {
    my $oid = $ps[7];
    exists $h->{$oid} || die "no mapping for $oid\n";
    $ps[7] = $h->{$oid};
  }
  if($opt eq "tgt") {
    my $oid = $ps[2];
    exists $h->{$oid} || die "no mapping for $oid\n";
    $ps[2] = $h->{$oid};
  }
  print $fho join(" ", @ps)."\n";
}
close $fhi;
close $fho;

__END__
