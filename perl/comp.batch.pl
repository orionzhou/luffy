#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.batch.pl - batch genome comparison

=head1 SYNOPSIS
  
  comp.batch.pl [-help] 

  Options:
    -h (--help)   brief help message

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
use Medicago;

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my ($tgt, @qrys) = ($tname, @qnames);
for my $qry (@qrys) {
  my $tag = $qry;
  $tag =~ s/HM//i;
  my $dir = "$ENV{'misc3'}/$qry\_$tgt/23_blat";
  chdir $dir || die "cannot chdir to $dir\n";
#  runCmd("comp.pl -q $qry -t $tgt");
#  runCmd("comp.novseq.pl -q $qry -t $tgt");
  runCmd("gal2chain.pl -i $dir/31.9.gal -o $dir/31.9.chain");
#  runCmd("comp.syn.ortho.pl -q $qry -t $tgt");
#  runCmd("comp.sv.pl -q $qry");
}
__END__

