#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.ortho.pl - 

=head1 SYNOPSIS
  
  ortho.group.pl [-help] assign proteins to ortholog groups

  Options:
    -h (--help)   brief help message
    -i (--in)     input fasta file
    -o (--out)    output (tabular) file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use Time::HiRes qw/gettimeofday tv_interval/;
use Tabix;
use Bio::DB::Fasta;
use Data::Dumper;
use Common;
use Location;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my $help_flag;
my ($fi, $fo) = ('', '');

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

#runCmd("qsub.blat.pl -n 1 -i $fi -o $fo.1 -t $fi -p pro");

### seq 0 23 | xargs -i printf "%02d\n" {} | parallel --no-notice -j 24 blat -prot $fi $fo.1/part.{}.fas $fo.1/part.{}.psl

#runCmd("cat $fo.1/part.*.psl > $fo.2.psl");
#runCmd("psl2gal.pl -i $fo.2.psl | gal.filter.ortho.pl | gal2mcl.pl -o $fo.4.tbl");
#runCmd("\$soft/mcl/bin/mcl $fo.4.tbl -te 4 -I 5.0 --abc -o $fo.5.mcl");
parse_mcl("$fo.5.mcl", $fo);

sub parse_mcl {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/grp org gid/)."\n";
  my $grp = 1;
  while(<$fhi>) {
    chomp;
    my @ps = split "\t";
    for (@ps) {
      my ($org, $id);
      if(/\-/) {
        ($org, $id) = split /\-/;
      } else {
        ($org, $id) = ("", $_);
      }
      print $fho join("\t", $grp, $org, $id)."\n";
    }
    $grp ++;
  }
  close $fhi;
  close $fho;
}
__END__

