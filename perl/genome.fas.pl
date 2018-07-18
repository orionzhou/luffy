#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  genome.fas.pl - process genome fasta file

=head1 SYNOPSIS
  
  genome.fas.pl [-help] [-org organism]

  Options:
    -h (--help)   brief help message
    -g (--org)    genmome ID (organism) to process

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;

my ($org) = ('');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "org|g=s"  => \$org,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$org;

my $dir = "$ENV{'genome'}/$org";
chdir $dir || die "cannot chdir to $dir\n";

runCmd("seq.check.pl -i raw.fas -o raw.fix.fas");
-s "raw.fix.fas.index" && runCmd("rm raw.fix.fas.index"); 
-l "11_genome.fas" && runCmd("rm 11_genome.fas"); 
-s "11_genome.fas.index" && runCmd("rm 11_genome.fas.index"); 
runCmd("seq.rename.pl -i raw.fix.fas -p scf -o 11_genome.fas");

-s "ctg.raw.fas" && runCmd("seq.check.pl -i ctg.raw.fas -o ctg.fas");

runCmd("seqlen.py 11_genome.fas 15.sizes");
runCmd("awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {print \$1, 0, \$2}' 15.sizes > 15.bed");

runCmd("seqgap.pl -i 11_genome.fas -o 16.gap.bed -m 10");
##runCmd("bedToBigBed 16.gap.bed 15.sizes 16.gap.bb");
runCmd("bgzip -c 16.gap.bed > 16.gap.bed.gz");
runCmd("tabix -p bed 16.gap.bed.gz");


__END__

