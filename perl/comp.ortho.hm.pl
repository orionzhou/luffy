#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.ortho.pl - 

=head1 SYNOPSIS
  
  comp.ortho.pl [-help]

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
use File::Path qw/make_path remove_tree/;
use File::Basename;
use Time::HiRes qw/gettimeofday tv_interval/;
use Bio::DB::Fasta;
use Data::Dumper;
use Common;
use Location;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $dir = "$ENV{'misc3'}/comp.ortho.hm";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $tgt = "HM101";
my @qrys = qw/
  HM058 HM125 HM056 HM129 HM060
  HM095 HM185 HM034 HM004 HM050 
  HM023 HM010 HM022 HM324 HM340
/;
my @orgs = ($tgt, @qrys);

ortho_aln_prep("21.ids.tbl", "25.aln.cmd", "25_seq", "25_aln");

sub ortho_aln_prep {
  my ($fi, $fo, $ds, $da) = @_;
  -d $ds || make_path($ds);
  -d $da || make_path($da);
  my $ti = readTable(-in => $fi, -header => 1);
  my @cnames = $ti->header;
  my @orgs = @cnames[1..$#cnames];

  my $h;
  for my $org (@orgs) {
    my $fs = "$ENV{'genome'}/$org/51.fas";
    $h->{$org} = Bio::DB::Fasta->new($fs);
  }

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  for my $i (0..$ti->lastRow) {
    my ($idx, @gids) = $ti->row($i);
    my (@nids, @seqs);
    for my $j (0..$#gids) {
      my $gid = $gids[$j];
      ($gid eq "" || $gid eq "-") && next;
      my $org = $orgs[$j];
      push @nids, "$org-$gid";
      my $seq = $h->{$org}->seq($gid);
      $seq =~ s/X//ig;
      push @seqs, $seq;
    }
    if(@nids > 1) {
      open(my $fhs, ">$ds/$idx.fas") or die "cannot write $idx\n";
      print $fhs join("\n", map {">".$nids[$_]."\n".$seqs[$_]} 0..$#nids);
      close $fhs;

      print $fho "clustalo -i $ds/$idx.fas -o $da/$idx.fas --outfmt=fasta --force --full --full-iter\n";
    }
  }
  close $fho;
}

__END__

