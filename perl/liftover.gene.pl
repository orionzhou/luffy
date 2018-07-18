#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  liftover.gene.pl - lift over a gene from one assembly to another

=head1 SYNOPSIS
  
  liftover.gene.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gtb)
    -o (--out)    output file
    -r (--ref)    target ref-seq fasta
    -g (--gtb)    target annotation (*.gtb.gz)
    -x (--gax)    aln block file (*.gax.gz)
    -s (--snp)    aln snp file (*.snp.gz)
    -d (--idm)    aln idm file (*.idm.gz)

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
use Bio::DB::Fasta;
use Common;
use Location;
use Seq;
use Gal;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my ($fr, $fg, $fx, $fs, $fd) = ('') x 5;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "ref|r=s"  => \$fr,
  "gtb|g=s"  => \$fg,
  "gax|x=s"  => \$fx,
  "snp|s=s"  => \$fs,
  "idm|d=s"  => \$fd,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;
pod2usage(2) if !$fr || !$fx || !$fs || !$fd;

my ($fhi, $fho);
if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $db = Bio::DB::Fasta->new($fr);
my $gax = Tabix->new(-data => $fx);
my $snp = Tabix->new(-data => $fs);
#my $idm = Tabix->new(-data => $fd);

my ($cnt, $cntu) = (0, 0);
print $fho join("\t", qw/id chr beg end srd cat3 clen lent qid qpos/)."\n";
while(<$fhi>) {
  chomp;
  /^(id)|(\#)/ && next;
  my ($id, $par, $chr, $beg, $end, $srd, 
    $elocs, $ilocs, $clocs, $flocs, $tlocs, $phase, 
    $src, $conf, $cat1, $cat2, $cat3, $note) = split "\t";
  next if $cat1 ne "mRNA";
  die "no CDS-loc for $id\n" unless $clocs;
  my $rcloc = locStr2Ary($clocs); 
  my $riloc = locStr2Ary($ilocs);
  my $cloc = $srd eq "+" ? 
    [ map {[$beg + $_->[0] - 1, $beg + $_->[1] - 1]} @$rcloc ] :
    [ map {[$end - $_->[1] + 1, $end - $_->[0] + 1]} @$rcloc ]; 
  my $iloc = $srd eq "+" ? 
    [ map {[$beg + $_->[0] - 1, $beg + $_->[1] - 1]} @$riloc ] :
    [ map {[$end - $_->[1] + 1, $end - $_->[0] + 1]} @$riloc ];
  my $staloc = $srd eq "+" ?
    [ $beg + $rcloc->[0]->[0] - 1, $beg + $rcloc->[0]->[0] + 1 ] : 
    [ $end - $rcloc->[0]->[0] - 1, $end - $rcloc->[0]->[0] + 1 ];
  my $stoloc = $srd eq "+" ?
    [ $beg + $rcloc->[-1]->[1] - 3, $beg + $rcloc->[-1]->[1] - 1 ] : 
    [ $end - $rcloc->[-1]->[1] + 1, $end - $rcloc->[-1]->[1] + 3 ];
#  my $seqsta = seqret_simple($db, $chr, @$staloc, $srd);
#  my $seqsto = seqret_simple($db, $chr, @$stoloc, $srd);
  
  my $clen = locAryLen($cloc);

  my $hs = {};
  for (@$cloc) {
    my ($ib, $ie) = @$_;
    my $ilen = $ie - $ib + 1;
    
    my $ary = read_gax($gax, $chr, $ib, $ie);
    for (@$ary) {
      my ($cid, $tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd, $lev) = @$_;
      $hs->{$cid} ||= [$lev, $qid, [], [], 0, 0];
      push @{$hs->{$cid}->[2]}, [$qb, $qe];
      $hs->{$cid}->[4] += $qe - $qb + 1;
    }
  }

  my $lent = scalar(keys(%$hs)) > 0 ? sum(map {$_->[4]} values(%$hs)) : 0;
  my @cids = sort {$hs->{$b}->[4] <=> $hs->{$a}->[4]} keys(%$hs);

  my ($qid, $qpos);
  if($lent / $clen >= 0.5) {
    my $cid = $cids[0];
    $qid = $hs->{$cid}->[1];
    my $qloc = $hs->{$cid}->[2];
    $qloc = [ sort {$a->[0] <=> $b->[0]} @$qloc ];
    $qpos = int( ($qloc->[0]->[0] + $qloc->[-1]->[1]) / 2 );
  } else {
    $cntu ++;
    $qid = "chrZ",
    $qpos = $cntu * 200000 + 1;
  }
  print $fho join("\t", $id, $chr, $beg, $end, $srd, $cat3,
    $clen, $lent, $qid, $qpos)."\n";
  $cnt ++;
}
close $fhi;
close $fho;
printf "%d out of %d un-aligned\n", $cntu, $cnt;

__END__
$ary = read_idm($idm, $chr, $ib, $ie);
    for (@$ary) {
      my ($tid, $tb, $te, $qid, $qb, $qe, $cid, $lev) = @$_;
      exists $hs->{$cid} || next;
      my ($ref, $olen) = posOvlp([[$tb+1, $te-1]], \@locs);
      $olen == 0 || next;
      
      $tb < $ie || next;
      $te > $ib || next;
      $tb = $ib if $tb < $ib;
      $te = $ie if $te > $ie;
      push @{$hs->{$cid}->[3]}, ($tb, $te);
    }

