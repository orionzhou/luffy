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
use Tabix;
use Bio::DB::Fasta;
use Data::Dumper;
use Common;
use Location;
use Gtb;
use Gal;
use Comp;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $dir = "$ENV{'misc3'}/comp.og";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $tgt = "HM101";
my @qrys = qw/
  HM058 HM125 HM056 HM129 HM060
  HM095 HM185 HM034 HM004 HM050 
  HM023 HM010 HM022 HM324 HM340
/;
my @orgs = ($tgt, @qrys);

#get_seq_by_org("01.gid.tbl", "02.fas");
### run comp.og.R

#get_seq_by_org("11.novel.gid.tbl", "12.fas");
#runCmd("qsub.blat.pl -n 1 -i 12.fas -o 15_blat -t $dir/12.fas -p pro");

### seq 0 23 | xargs -i printf "%02d\n" {} | parallel --no-notice -j 24 blat -prot 12.fas 15_blat/part.{}.fas 15_blat/part.{}.psl

#runCmd("cat 15_blat/part.*.psl > 15.blat.psl");
#runCmd("psl2gal.pl -i 15.blat.psl | gal.filter.ortho.pl | gal2mcl.pl -o 16.tbl");
#runCmd("\$soft/mcl/bin/mcl 16.tbl -te 4 -I 5.0 --abc -o 16.mcl");
#parse_mcl("16.mcl", "17.group.tbl");

sub get_seq_by_org {
  my ($fi, $fo) = @_;
  my $hi;
  my $ti = readTable(-in => $fi, -header => 1);
  for my $i (0..$ti->lastRow) {
    my ($org, $gid, $idx) = $ti->row($i);
    $hi->{$org} ||= [];
    push @{$hi->{$org}}, [$gid, $idx];
  }

  my $seqHO = Bio::SeqIO->new(-file => ">$fo", -format => 'fasta');
  my @orgs = sort(keys(%$hi));
  for my $org (@orgs) {
    my $fs = "$ENV{'genome'}/$org/51.fas";
    my $db = Bio::DB::Fasta->new($fs);
    for (@{$hi->{$org}}) {
      my ($gid, $idx) = @$_;
      my $seq = $db->seq($gid);
      my $nid = "$org-$gid";
      #my $nid = defined($idx) ? "$org-$gid-$idx" : "$org-$gid";
      $seqHO->write_seq(Bio::Seq->new(-id => $nid, -seq => $seq));
    }
  }
  $seqHO->close();
}
sub write_one_hit {
  my ($tid, $h, $fho) = @_;
  for my $org (keys(%$h)) {
    my ($qid, $score) = @{$h->{$org}};
    print $fho join("\t", $tid, $qid, $score)."\n";
  }
}
sub get_best_hit {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  my $h;
  my $ptid = "";
  while(<$fhi>) {
    chomp;
    my ($tid, $qid, $torg, $qorg, $e1, $e2, $ident, $cov) = split "\t";
    next if $torg eq $qorg || $cov < 0.5;
    my $score = $ident * $cov / 10000;
    if($ptid ne $tid && $ptid ne "") {
      write_one_hit($ptid, $h, $fho);
      $h = {$qorg => [$qid, $score]};
    } else {
      $h->{$qorg} ||= [$qid, $score];
      $h->{$qorg} = [$qid, $score] if $score > $h->{$qorg}->[1];
    }
    $ptid = $tid;
  }
  write_one_hit($ptid, $h, $fho);
  close $fhi;
  close $fho;
}
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
      my ($org, $id) = split /\-/;
      print $fho join("\t", $grp, $org, $id)."\n";
    }
    $grp ++;
  }
  close $fhi;
  close $fho;
}
__END__

