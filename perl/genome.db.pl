#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  genome.db.pl - make genome db

=head1 SYNOPSIS
  
  genome.db.pl [-help] [-org organism]

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

makedb_blat();
#makedb_blast();
#makedb_bwa();
makedb_bowtie2();

sub makedb_blat {
  my $dir = "$ENV{'data'}/db/blat";
  chdir $dir || die "cannot chdir to $dir\n";
  runCmd("faToTwoBit $fg $org.2bit");
  runCmd("blat $org.2bit tmp.fas tmp.out -makeOoc=$org.2bit.tile11.ooc");
  -s "tmp.out" && runCmd("rm tmp.out");
}
sub makedb_blast {
  my $dir = "$ENV{'data'}/db/blast";
  chdir $dir || die "cannot chdir to $dir\n";
  runCmd("makeblastdb -dbtype nucl -in $fg -out $org -title $org");
}
sub makedb_bwa {
  my $dir = "$ENV{'data'}/db/bwa";
  chdir $dir || die "cannot chdir to $dir\n";
  runCmd("bwa index -a is -p $dir/$org $fg");
}
sub makedb_bowtie2 {
  my $dir = "$ENV{'data'}/db/bowtie2";
  chdir $dir || die "cannot chdir to $dir\n";
  -s "$org.fas" && runCmd("rm $org.fas");
  runCmd("ln -sf $fg $org.fa");
  runCmd("bowtie2-build $org.fa $org");
}
__END__

