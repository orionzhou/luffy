#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  update.taxdb.pl - update local NCBI Taxonomy DB

=head1 SYNOPSIS
  
  update.taxdb.pl [-help]

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
use Common;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $dir = "$ENV{'data'}/db/ncbi_taxon";
chdir $dir || die "cannot chdir to $dir\n";

my $url_pre = "ftp://ftp.ncbi.nih.gov//pub/taxonomy";
my @fns;

@fns = qw/gi_taxid_nucl.dmp.gz gi_taxid_prot.dmp.gz/;
for my $fn (@fns) {
#  -s $fn && runCmd("rm $fn");
#  runCmd("wget $url_pre/$fn");
#  runCmd("gunzip $fn");
}

@fns = qw/taxdump.tar.gz/;
for my $fn (@fns) {
  -s $fn && runCmd("rm $fn");
  runCmd("wget $url_pre/$fn");
  runCmd("tar xzf $fn");
}


exit 0;
