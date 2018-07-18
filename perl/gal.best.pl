#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal.best.pl - pick best hits for query/target

=head1 SYNOPSIS
  
  gal.best.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -p (--opt)    option (tgt | qry, def: qry)
    -r (--rep)    whether discard repetitive mapping (def: FALSE)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use Gal;

my ($fi, $fo) = ('') x 2;
my $opt = 'qry';
my $flag_rep = 0;
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "opt|p=s"  => \$opt,
  "rep|r"    => \$flag_rep,
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

my ($hq, $ht);
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps == 21;
  my ($cid, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $lev, $ali, $mat, $mis, $qN, $tN, $ident, $score, $tlS, $qlS) = @$ps;
  my ($qLoc, $tLoc) = (locStr2Ary($qlS), locStr2Ary($tlS));
  @$qLoc == @$tLoc || die "unequal pieces\n";
  my $nBlock = @$qLoc;

  $hq->{$qId} ||= [];
  $ht->{$tId} ||= [];
  push @{$hq->{$qId}}, [$ps, $score];
  push @{$ht->{$tId}}, [$ps, $score];
}
close $fhi;

die "unknown opt: $opt\n" if $opt ne "qry" && $opt ne "tgt";
my $h = $opt eq "qry" ? $hq : $ht;

print $fho join("\t", @HEAD_GAL)."\n";
for my $uid (sort(keys(%$h))) {
  my @ps = @{$h->{$uid}};
  @ps = sort {$b->[1] <=> $a->[1]} @ps;
  next if $flag_rep && @ps >= 2 && $ps[0]->[1] == $ps[1]->[1];
  print $fho join("\t", @{$ps[0]->[0]})."\n";
}
close $fho;


__END__
