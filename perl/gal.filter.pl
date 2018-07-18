#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal.filter.pl - filter Gal records 

=head1 SYNOPSIS
  
  gal.filter.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -m (--match)  minimum matches (default: 1) 
    -p (--ident)  minimum percent identity (default: 0.5) 
    -s (--score)  minimum score (default: 0) 
    -c (--cov)    minimum query coverage (default: 0.1) 

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
my ($min_match, $min_ident, $min_score, $min_qcov) = (1, 0.5, 0, 0);
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "match|m=i"  => \$min_match,
  "ident|p=f"  => \$min_ident,
  "score|s=i"  => \$min_score,
  "cov|c=f"    => \$min_qcov,
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

print $fho join("\t", @HEAD_GAL)."\n";

my $cnt = 0;
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [split "\t"];
  next unless @$ps == 21;
  my ($cid, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $lev, $ali, $mat, $mis, $qN, $tN, $ident, $score, $tlS, $qlS) = @$ps;
  $mat >= $min_match || next;
  $ident >= $min_ident || next;
  $ali / $qSize >= $min_qcov || next;
  next if $score ne "" && $score < $min_score;
  print $fho join("\t", @$ps)."\n";
  $cnt ++;
}
print STDERR "$cnt rows passed filter\n";
close $fhi;
close $fho;


__END__
