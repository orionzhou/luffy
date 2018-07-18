#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use Getopt::Long;
use Pod::Usage;
use Common;
use Data::Dumper;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my ($fhi, $fho);
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "Can't open file $fi for reading: $!\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $hi;
print $fho "##gff-version 3\n";
while(<$fhi>) {
  chomp;
  next if !$_ || /^\#/;
  my @ps = split "\t";

#  my $h = parse_gff_tags($ps[8]);
#  if(exists $h->{"ID"}) {
#    $h->{"ID"} = change_augus_id($h->{"ID"});
#  }
#  if(exists $h->{"Parent"}) {
#    $h->{"Parent"} = change_augus_id($h->{"Parent"});
#  }
#  $ps[8] = join(";", map {$_."=".$h->{$_}} sort(keys(%$h)));
  
  $ps[2] =~ /^(gene)|(transcript)|(CDS)$/ || next;
  if($ps[2] eq "transcript") {
    $ps[2] = "mRNA";
    $ps[8] .= ";Note=$ps[5]";
  }
  $ps[5] = '.';
  print $fho join("\t", @ps)."\n";
}
close $fhi;
close $fho;

sub change_augus_id {
  my ($oid) = @_;
  my $nid = 'err';
  if($oid =~ /^g(\d+)$/) {
    $nid = sprintf "g%05d", $1;
  } elsif($oid =~ /^g(\d+)\.t(\d+)$/) {
    $nid = sprintf "g%05d.%d", $1, $2;
  } elsif($oid =~ /^g(\d+)\.t(\d+)\.(\w+)/) {
    $nid = sprintf "g%05d.%d.%s", $1, $2, $3;
  } else {
    die "unsupported ID: $oid\n";
  }
  return $nid;
} 
