#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  genome.anno.pl - 

=head1 SYNOPSIS
  
  genome.anno.pl [-help] [-org organism]

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

my $fg = "$ENV{'genome'}/$org/11_genome.fas";
-s $fg || die "$fg is not there\n";

### needs to run mt.rnaseq.pl on Mesabi first
### needs to be run on Mesabi
#runCmd("repeatmasker.pl -g $org");
#runCmd("mt.augus.pl -g $org");

runCmd("mt.nbs.pl -g $org");
runCmd("mt.rlk.pl -g $org");
runCmd("spada.pl --cfg \$spada/cfg.txt \\
  --dir \$misc4/spada.crp.$org \\
  --hmm \$misc4/hmm/crp \\
  --fas \$genome/$org/11_genome.fas \\
  --gff \$genome/$org/augustus/31.gff \\
  --org Mtruncatula --sp --threads 24");
runCmd("mt.anno.pl -g $org");

__END__

