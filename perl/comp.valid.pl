#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.valid.pl - 

=head1 SYNOPSIS
  
  comp.valid.pl [-help] 

  Options:
    -h (--help)   brief help message

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Common;
use Location;
use List::Util qw/min max sum/;

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $dir = "$ENV{'misc3'}/seqvalidation";
chdir $dir or die "cannot chdir to $dir\n";

my $qry_fas = "$dir/01.fas";
my $tgt_fas = "$ENV{'genome'}/HM101/11_genome.fas";

#runCmd("seq.check.pl -i 00.raw.fas -o 01.fas");
#runCmd("qsub.blat.pl -i 01.fas -o 05");
## submit itasca
#runCmd("cat 05/part.*.psl > 11.psl");
#runCmd("psl2gal.pl -i 11.psl -o 11.gal");
#runCmd("gal.fix.ovlp.pl -i 11.gal -o - | \\
#  gal.calib.pl -i - -q $qry_fas -t $tgt_fas -o 12.fixed.gal");
#runCmd("gal.best.pl -r -i 12.fixed.gal -o 15.best.gal");
runCmd("gal2snp.pl -i 15.best.gal -o 16.snp -q $qry_fas -t $tgt_fas");

__END__
