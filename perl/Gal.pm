package Gal;
use strict;
use Data::Dumper;
use Common;
use Location;
use Seq;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/@HEAD_GAL $h_score
    read_gax 
    read_snp read_idm read_idm_bp/;
@EXPORT_OK = qw//;

our @HEAD_GAL = qw/id tId tBeg tEnd tSrd tSize
  qId qBeg qEnd qSrd qSize
  lev ali mat mis qN tN ident score tLoc qLoc/;

our $h_score = { # lev -> score
  1 => 1000,
  2 => 500,
  3 => 250,
  4 => 125,
  5 => 63,
};

  my ($cid, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $lev, $ali, $mat, $mis, $qN, $tN, $ident, $score, $tlS, $qlS) = ();

sub read_gax {
  my ($con, $chr, $beg, $end) = @_;
  my $iter = $con->query($chr, $beg - 1, $end);
  my @ary;
  return \@ary if ! $iter->get();
  while (my $line = $con->read($iter)) {
    my @ps = split("\t", $line);
    my ($tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd, $cid, $lev) = @ps;
    $te - $tb == $qe - $qb || die "not gax: $tid:$tb-$te != $qid:$qb-$qe\n";
    
    if($tb < $beg) {
      $qb += $beg - $tb if $tsrd eq $qsrd;
      $qe -= $beg - $tb if $tsrd ne $qsrd;
      $tb = $beg;
    } 
    if($te > $end) {
      $qe -= $te - $end if $tsrd eq $qsrd;
      $qb += $te - $end if $tsrd ne $qsrd;
      $te = $end;
    }
#    my $alen = $te - $tb + 1;
#    print "$tid : $tb - $te, $qid : $qb - $te, $alen\n";
    push @ary, [$cid, $tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd, $lev];
  }
  return \@ary;
}
sub read_snp {
  my ($con, $id, $beg, $end) = @_;
  my $iter = $con->query($id, $beg - 1, $end);
  my @snps;
  while (my $line = $con->read($iter)) {
    my ($tid, $tpos, $tnt, $qnt, $qid, $qpos, $cid, $lev) = 
      split("\t", $line);
    next if $tpos < $beg || $tpos > $end;
    push @snps, [$tid, $tpos];
  }
  return \@snps;
}
sub read_idm {
  my ($con, $id, $beg, $end) = @_;
  my $iter = $con->query($id, $beg - 1, $end);
  my @idms;
  while (my $line = $con->read($iter)) {
    my ($tid, $tb, $te, $qid, $qb, $qe, $cid, $lev) = split("\t", $line);
    push @idms, [$tid, $tb, $te, $qid, $qb, $qe, $cid, $lev];
  }
  return \@idms;
}
sub read_idm_bp {
  my ($con, $id, $beg, $end, $max_len) = @_;
  my $iter = $con->query($id, $beg - 1, $end);
  my $h = {};
  while (my $line = $con->read($iter)) {
    my ($tid, $tb, $te, $qid, $qb, $qe, $cid, $lev) = split("\t", $line);
    $te - $tb - 1 <= $max_len || next;
    $qe - $qb - 1 <= $max_len || next;
    $tb < $end || next;
    $te > $beg || next;
    $tb = $beg if $tb < $beg;
    $te = $end if $te > $end;
    my ($del, $ins) = ($te - $tb - 1, $qe - $qb - 1);
    $h->{$cid} ||= [0, 0, $lev];
    $h->{$cid}->[0] += $del;
    $h->{$cid}->[1] += $ins;
  }
  return $h;
}

1;
__END__
