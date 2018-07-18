#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gax2bed.pl - convert a Gax (long) file to BED file

=head1 SYNOPSIS
  
  gax2bed.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -p (--opt)    option ('qry' or 'tgt', default: 'tgt')

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gal;

my ($fi, $fo) = ('') x 2;
my $opt = "tgt";
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "opt|p=s"  => \$opt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$opt;
$opt = lc($opt);

if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot wriet $fo\n";
}

while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my ($tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd, $cid, $lev) 
    = split "\t";
  my $srd = is_revsrd($tsrd, $qsrd) ? "-" : "+";
  my $score = exists $h_score->{$lev} ? $h_score->{$lev} : 300;
  if($opt eq "qry") {
    my $id = "$cid|$tid:$tb-$te";
    print $fho join("\t", $qid, $qb - 1, $qe, $id, $score, $srd)."\n";
  } elsif($opt eq "tgt") {
    my $id = "$cid|$qid:$qb-$qe";
    print $fho join("\t", $tid, $tb - 1, $te, $id, $score, $srd)."\n";
  } else {
    die "unknow optiotn: $opt\n";
  }
}
close $fhi;
close $fho;


__END__
