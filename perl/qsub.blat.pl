#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  qsub.blat.pl - split a fasta file and build qsub commands

=head1 SYNOPSIS
  
  qsub.blat.pl [-help] [-in input-file] [-out output-directory]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (fasta) file
    -o (--out)    output directory
    -n (--num)    number of qsub batches (def: 1)
    -t (--tgt)    blat target (def: "HM101")
    -g (--tag)    qsub job tag (def: "pz")
    -p (--opt)    which qsub script to run (1: blat; 2: blat2; 3: blat-pro)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use File::Path qw/make_path remove_tree/;
use File::Spec;
use Common;
use Data::Dumper;

my ($fi, $dir) = ('') x 2;
my ($n, $tgt) = (1, "HM101");
my ($tag, $opt) = ("pz", 1);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$dir,
  "num|n=i" => \$n,
  "tgt|t=s" => \$tgt,
  "tag|g=s" => \$tag,
  "opt|p=s" => \$opt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$dir;

$fi = File::Spec->rel2abs($fi);
$dir = File::Spec->rel2abs($dir);
$tgt = File::Spec->rel2abs($tgt) if $opt eq "pro";

runCmd("rm -rf $dir") if -d $dir;
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

runCmd("ln -sf $fi part.fas");

my $ppn = 24;
my $part = $n * $ppn;
my $digits = getDigits($part);
runCmd("pyfasta split -n $part part.fas");
runCmd("rm part.fas part.fas.*");

my $ps = runCmd("du -h part.*.fas", 2);
my @sizes;
for (@$ps) { push @sizes, [ split " ", $_ ]->[0]; }
@sizes = sort @sizes;
print "\n##### stats begin #####\n";
printf "range: %s  -  %s\n", $sizes[0], $sizes[$#sizes];
print "##### stats end   #####\n\n";

print "\n##### qsub command begins #####\n";
my $prog;
if($opt == 1) {
  $prog = "comp2";
} elsif($opt == 2) {
  $prog = "blat2";
} elsif($opt == 3) {
  $prog = "blat-pro";
} else {
  die "unknonw opt: $opt\n";
}

if($n == 1) {
  print "qsub $prog -N blat.$tag -t 1 -v PRE=$dir/part,SUF=fas,DIG=$digits,TGT=$tgt\n";
} else {
  print "qsub $prog -N blat.$tag -t 1-$n -v PRE=$dir/part,SUF=fas,DIG=$digits,TGT=$tgt\n";# -l qos=weightlessqos\n";
}
print "##### qsub command ends   #####\n\n";



exit 0;
