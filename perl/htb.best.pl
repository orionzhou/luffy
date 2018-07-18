#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  htb.best.pl - pick best hmm in Htb

=head1 SYNOPSIS
  
  htb.best.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use Htb;

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
#pod2usage(2) if !$fi || !$fo;

if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $h = {};
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [split("\t", $_, -1)];
  next unless @$ps == 15;
  my ($id, $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $hId, $hBeg, $hEnd, $hSrd, $hSize, 
    $e, $score, $qlS, $hls) = @$ps;
  if(exists($h->{$hId})) {
    my $pscore = $h->{$hId}->[12];
    $h->{$hId} = $ps if $pscore < $score;
  } else {
    $h->{$hId} = $ps;
  }
}
close $fhi;

print $fho join("\t", @HEAD_HTB)."\n";
for my $hId (sort(keys(%$h))) {
  my $ps = $h->{$hId};
  print $fho join("\t", @$ps)."\n";
}
print STDERR scalar(keys(%$h))." rows picked\n";
close $fho;


__END__
