#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  genomestat.pl - genome statistics

=head1 SYNOPSIS
  
  genomestat.pl [-help] [-o organism]

  Options:
    -h (-help)   brief help message
    -o (--org)   organism

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw/make_path remove_tree/;
use Common;

my ($org) = ('');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "org|o=s"  => \$org,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$org;

my $dir = "/home/youngn/zhoup/Data/genome/$org";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

wc -l 51.gtb
mergeBed -i 51.bed/mrna.bed | bedlen.pl

wc -l 43_crp.gtb
mergeBed -i 43_crp.bed/mrna.bed | bedlen.pl

wc -l 42_nbs.gtb
mergeBed -i 42_nbs.bed/mrna.bed | bedlen.pl

__END__
