#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal.stat.pl -

=head1 SYNOPSIS
  
  gal.stat.pl [-help] [-in input-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gal)

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
use Gal;

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
pod2usage(2) if !$fi;

runCmd("awk 'BEGIN{OFS=\"\\t\"} {if(NR > 1) print \$2, \$3-1, \$4}' $fi > $fi.t.bed");
runCmd("awk 'BEGIN{OFS=\"\\t\"} {if(NR > 1) print \$7, \$8-1, \$9}' $fi > $fi.q.bed");
runCmd("sortBed -i $fi.t.bed | mergeBed -i stdin > $fi.t.r.bed");
runCmd("sortBed -i $fi.q.bed | mergeBed -i stdin > $fi.q.r.bed");

my $tids = runCmd("awk '{if(NR>1 && NR<11) print \$2}' $fi | sort | uniq", 2);
my $qids = runCmd("awk '{if(NR>1 && NR<11) print \$7}' $fi | sort | uniq", 2);
my $tlena = runCmd("bedlen.pl -i $fi.t.bed", 2)->[0];
my $qlena = runCmd("bedlen.pl -i $fi.q.bed", 2)->[0];
my $tlenr = runCmd("bedlen.pl -i $fi.t.r.bed", 2)->[0];
my $qlenr = runCmd("bedlen.pl -i $fi.q.r.bed", 2)->[0];

printf "tgt ids: %s\n", join(" ", @$tids);
printf "  total   len: %10d\n", $tlena;
printf "  reduced len: %10d\n", $tlenr;
printf "qry ids: %s\n", join(" ", @$qids);
printf "  total   len: %10d\n", $qlena;
printf "  reduced len: %10d\n", $qlenr;

runCmd("rm $fi.t.bed $fi.q.bed $fi.t.r.bed $fi.q.r.bed", -1);
__END__
