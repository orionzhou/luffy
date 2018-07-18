#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  usearch.pl - cluster sequences using USEARCH

=head1 SYNOPSIS
  
  novseq.pl [-help] [-in input-fasta] [-o output-prefix]

  Options:
    -h (--help)   brief help message
    -i (--in )    input sequence file (fasta)
    -o (--out)    output prefix
    -r (--red)    percent redundancy threshold (0.001, i.e., 0.1%)
    -l (--len)    minimum length of a cluster to keep (50)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Data::Dumper;
use Common;
use Location;
use Align;
use Seq;
use File::Basename;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fi, $fo) = ('') x 2;
my $red = 0.001;
my $minlen = 50;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "red|r=f" => \$red,
  "len|l=i" => \$minlen,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my $mincols = 50;
my ($redcy, $rep) = (1, 0);

my $dir = $fo;
-d $dir || make_path($dir);

runCmd("usearch -sortbylength $fi -output $dir/r0.fas");
runCmd("seqlen.pl -i $dir/r0.fas | awk 'BEGIN{OFS=\"\\t\"; C=1} {print \$1, 1, \$2, \"+\", C; C+=1}' > $dir/r0.clu");

while($redcy > $red) {
  printf "=============== Begin Rep %d ===============\n", ++$rep;
#  next if $rep <= 1;
  my $prev_fas = sprintf("%s/r%d.fas", $dir, $rep - 1);
  my $prev_clu = sprintf("%s/r%d.clu", $dir, $rep - 1);
  my $pre = sprintf("%s/r%d", $dir, $rep);
  runCmd("usearch -cluster_smallmem $prev_fas \\
    -id 0.9 -mincols $mincols -strand both -uc $pre.1.uc");
  uc_parse("$pre.1.uc", "$pre.2.tbl", $mincols, $rep);
  uc_split("$pre.2.tbl", "$pre.3.clu");
  my $nbp_red = cluster_stat("$pre.3.clu");
  my $nbp = get_fas_len($prev_fas);
  $redcy = $nbp_red / $nbp;

  update_cluster($prev_clu, "$pre.clu", "$pre.3.clu");
  cluster_stat("$pre.clu");
  clu2aln("$pre.clu", $fi, "$pre.5.aln");
  clu2bed("$pre.clu", "$pre.bed", $minlen);
  runCmd("seqret.pl -d $fi -b $pre.bed -o $pre.fas");
  runCmd("rm $pre.2.tbl $pre.3.clu");
  printf "##### redundancy: %.05f\n", $redcy;
}
runCmd("ln -sf $dir/r$rep.fas $fo.fas");
runCmd("ln -sf $dir/r$rep.bed $fo.bed");
runCmd("ln -sf $dir/r$rep.clu $fo.clu");

sub read_cluster {
  my ($fi) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  my ($hl, $hc);
  while(<$fhi>) {
    chomp;
    my ($id, $beg, $end, $srd, $cid) = split "\t";
    $hl->{$id} ||= [];
    push @{$hl->{$id}}, [$beg, $end, $srd, $cid];
    my $idx = @{$hl->{$id}} - 1;
    $hc->{$cid} ||= [];
    push @{$hc->{$cid}}, [$id, $idx];
  }
  close $fhi;
  return ($hl, $hc);
}
sub write_cluster {
  my ($hl, $fo) = @_;
  open(my $fho, ">$fo") or die "cannot write to $fo\n";
  for my $sid (sort(keys(%$hl))) {
    my @locs = @{$hl->{$sid}};
    for (@locs) {
      my ($beg, $end, $srd, $cid) = @$_;
      print $fho join("\t", $sid, $beg, $end, $srd, $cid)."\n";
    }
  }
  close $fho;
}
sub cluster_stat {
  my ($fi) = @_;
  my ($hl, $hc) = read_cluster($fi);
  
  my ($nc, $ns, $ne) = (0, 0, 0);
  my ($lt, $li, $ln, $lc) = (0, 0, 0, 0);
  for my $sid (keys(%$hl)) {
    $ns ++;
    my @locs = @{$hl->{$sid}};
    $ne += @locs;
  }
  for my $cid (keys(%$hc)) {
    my @refs = @{$hc->{$cid}};
    my $len;
    $nc ++;
    for (@refs) {
      my ($id, $idx) = @$_;
      my ($beg, $end, $srd) = @{$hl->{$id}->[$idx]};
      $len ||= $end - $beg + 1;
      $len == $end - $beg + 1 || die "len err: $id: $beg - $end != $len\n";
    }
    $lt += @refs * $len;
    if(@refs == 1) {
      $li += $len;
    } else {
      $lc += $len * @refs;
      $ln += $len;
    }
  }

  printf "seqs [%6d] segs [%6d] clst [%6d]\n", $ns, $ne, $nc;
  printf "total [%10d]: sgnt [%10d] + clsr    [%10d]\n", $lt, $li, $lc;
  printf "dedup [%10d]: sgnt [%10d] + clsr-nr [%10d]\n", $li+$ln, $li, $ln;
  $lt == $li + $lc || die "len error\n";
  return $lc - $ln;
}
sub update_cluster {
  my ($fi, $fo, $fl) = @_;
  my ($hl, $hc) = read_cluster($fi);
  my ($hll, $hcl) = read_cluster($fl);

  for my $lcid (keys(%$hcl)) {
    my @refs = @{$hcl->{$lcid}};
    next if @refs == 1;
    my @clu;
    for (@refs) {
      my ($id, $idx) = @$_;
      my ($beg, $end, $srd) = @{$hll->{$id}->[$idx]};
      push @clu, [$id, $beg, $end, $srd];
    }
    ($hl, $hc) = insert_cluster($hl, $hc, \@clu);
  }
  write_cluster($hl, $fo);
}
sub find_loc_idx {
  my ($hl, $id, $beg, $end) = @_;
  exists $hl->{$id} || die "no sid: $id in hl\n";
  my @locs = @{$hl->{$id}};

  my @idxs = indexes {$_->[0] <= $beg && $_->[1] >= $end} @locs;
  @idxs == 1 || die "$id: $beg-$end not part of loc\n".Dumper(\@locs);
  return $idxs[0];
}
sub insert_cluster {
  my ($hl, $hc, $clu) = @_;
  my $cid = max(keys(%$hc)) + 1;
  
  my $ncid = $cid + 1;
  for (@$clu) {
    my ($xi, $xnb, $xne, $xns) = @$_;
    my $xidx = find_loc_idx($hl, $xi, $xnb, $xne);
    
    my ($xob, $xoe, $xos, $ocid) = @{$hl->{$xi}->[$xidx]};

    for (@{$hc->{$ocid}}) {
      my ($yi, $yidx) = @$_;
      my ($yob, $yoe, $yos, $ocid2) = @{$hl->{$yi}->[$yidx]};
      $ocid == $ocid2 || die "hc/hl error: $yi:$yob-$yoe\n";
      my $yns = $xos eq $yos ? $xns : get_revsrd($xns);
      my $ynb = coordTransform($xnb, [[$xob, $xoe]], $xos, [[$yob, $yoe]], $yos);
      my $yne = coordTransform($xne, [[$xob, $xoe]], $xos, [[$yob, $yoe]], $yos);
      ($ynb, $yne) = ($yne, $ynb) if $ynb > $yne;

      if($yob < $ynb) {
        push @{$hl->{$yi}}, [$yob, $ynb-1, $yos, $ncid];
        $hc->{$ncid} = [[$yi, @{$hl->{$yi}}-1]];
        $ncid ++;
      }
      if($yne < $yoe) {
        push @{$hl->{$yi}}, [$yne+1, $yoe, $yos, $ncid];
        $hc->{$ncid} = [[$yi, @{$hl->{$yi}}-1]];
        $ncid ++;
      }
      
      $hl->{$yi}->[$yidx] = [$ynb, $yne, $yns, $cid];
      $hc->{$cid} ||= [];
      push @{$hc->{$cid}}, [$yi, $yidx];
    }
    delete $hc->{$ocid};
  }
  return ($hl, $hc);
}
sub parse_cigar {
  my ($str) = @_;
  my (@tloc, @qloc);
  my ($tp, $qp) = (1, 1);
  while($str =~ /(\d*)([A-Z])/g) {
    my $len = $1 eq "" ? 1 : $1;
    if($2 eq "I") {
      $tp += $len;
    } elsif($2 eq "D") {
      $qp += $len;
    } elsif($2 eq "M") {
      push @tloc, [$tp, $tp + $len - 1];
      push @qloc, [$qp, $qp + $len - 1];
      $tp += $len;
      $qp += $len;
    } else {
      die "unknow cigar character: $2 in $str\n";
    }
  }
  return (\@tloc, \@qloc);
}
sub recover_cigar {
  my ($tid, $qid, $db) = @_;
  my @ops;
  my $tseq = $db->seq($tid);
  my $qseq = $db->seq($qid);
  my ($tlen, $qlen) = (length($tseq), length($qseq));
  if($tseq =~ /^(\w*)$qseq(\w*)$/) {
    push @ops, ["I", length($1)] if $1 ne "";
    push @ops, ["M", $qlen];
    push @ops, ["I", length($2)] if $2 ne "";
  } elsif($qseq =~ /^(\w*)$tseq(\w*)$/) {
    push @ops, ["D", length($1)] if $1 ne "";
    push @ops, ["M", $tlen];
    push @ops, ["D", length($2)] if $2 ne "";
  } else {
    die "$tid [$tseq]\n$qid [$qseq] not match\n";
  }
  return join("", map {$_->[1].$_->[0]} @ops);
}
sub change_coord {
  my ($id, $beg, $end) = @_;
  $id =~ /^([\w\-\.]+)\-([0-9e\+]+)\-([0-9e\+]+)$/ || die "err id: $id\n";
  my ($nid, $ob, $oe) = ($1, $2, $3);
  my $nb = $ob + $beg - 1;  
  my $ne = $ob + $end - 1;  
  return ($nid, $nb, $ne);
}
sub uc_parse {
  my ($fi, $fo, $mincols, $rep) = @_;
#  my $db = Bio::DB::Fasta->new($fs);

  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  my ($h, $hl, $hc);
  while(<$fhi>) {
    chomp;
    my ($type, $id, $len, $idt, $srd, $t1, $t2, $aln, $qid, $tid) = 
      split "\t";
    if($type eq "S") {
      $hl->{$qid} = $len;
      $h->{$qid} = [];
    } elsif($type eq "H") {
      exists $hl->{$tid} || die "$tid not found\n";
      push @{$h->{$tid}}, [$qid, $srd, $len, $idt, $aln];
      $hc->{$qid} = $aln;
    } elsif($type eq "C") {
      $len - 1 == @{$h->{$qid}} || die "not $len children: $qid\n";
    }
  }
  close $fhi;
  my $nc = scalar(keys(%$h));
  my $ns = sum( map {scalar(@{$_})+1} values(%$h) );
  print "\t$nc clus: $ns seqs\n";

  for my $tid (sort(keys(%$h))) {
    my $tlen = $hl->{$tid};
    my ($ntid, $ntb, $nte) = $rep == 1 ? 
      ($tid, 1, $tlen) : change_coord($tid, 1, $tlen);
    
    my @strs;
    for (@{$h->{$tid}}) {
      my ($qid, $qsrd, $qlen, $idt, $aln) = @$_;
      $aln !~ /^\=/ || die "error cigar: $qid: $aln\n";

      my ($tloc, $qloc) = parse_cigar($aln);
      my @idxs = indexes {$_->[1] - $_->[0] + 1 >= $mincols} @$tloc;
      next if @idxs == 0;

      $tloc = [ @$tloc[@idxs] ];
      $qloc = [ @$qloc[@idxs] ];
      $qloc = [ reverse map { [$qlen-$_->[1]+1, $qlen-$_->[0]+1] } 
        @$qloc ] if $qsrd eq "-";
      locAryLen($qloc) == locAryLen($tloc) || die "$tid loc error\n";
      my ($tls, $qls) = (locAry2Str($tloc), locAry2Str($qloc));
      
      my ($nqid, $nqb, $nqe) = $rep == 1 ? 
        ($qid, 1, $qlen) :  change_coord($qid, 1, $qlen);
      my $str = join(":", $nqid, $nqb, $nqe, $qsrd, $tls, $qls);
      push @strs, $str;
    }
    next if @strs == 0;
    my $tstr = join(":", $ntid, $ntb, $nte, "+");
    print $fho join("\t", $tstr, join(" ", @strs))."\n";
  }
  close $fho;
}
sub uc_split {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";

  my $cid = 1;
  while(<$fhi>) {
    chomp;
    my ($tstr, $str) = split "\t";
    my ($tid, $tbeg, $tend, $tsrd) = split(":", $tstr);
    my @ps = split(" ", $str);
    next if @ps == 0;

    my (@qstats, @locs);
    for my $i (0..$#ps) {
      my ($qid, $qbeg, $qend, $qsrd, $tls, $qls) = split(":", $ps[$i]);
      my ($tloc, $qloc) = (locStr2Ary($tls), locStr2Ary($qls));
      $tloc = [ map {[$tbeg+$_->[0]-1, $tbeg+$_->[1]-1]} @$tloc ];
      $qloc = [ map {[$qbeg+$_->[0]-1, $qbeg+$_->[1]-1]} @$qloc ];
      push @qstats, [$qid, $qbeg, $qend, $qsrd, $tloc, $qloc];

      for my $j (0..@$tloc-1) {
        my ($tb, $te) = @{$tloc->[$j]};
        my ($qb, $qe) = @{$qloc->[$j]};
        die Dumper($qid, \@locs) if $qe > $qend;
        push @locs, [$tb, $te, "$i-$j"];
      }
    }
    
    my $ref = posSplit(\@locs);
    for (@$ref) {
      my ($beg, $end, $idxs) = @$_;
      my $len = $end - $beg + 1;
      print $fho join("\t", $tid, $beg, $end, "+", $cid)."\n";
      
      for my $idxS (@$idxs) {
        my ($idx) = split("-", $idxS);
        my ($qid, $qbeg, $qend, $qsrd, $tl, $ql) = @{$qstats[$idx]};
        my $qb = coordTransform($beg, $tl, $tsrd, $ql, $qsrd);
        my $qe = coordTransform($end, $tl, $tsrd, $ql, $qsrd);
        ($qb, $qe) = ($qe, $qb) if $qb > $qe;
        $qe - $qb == $end - $beg || 
          die "error: $tid\[$beg-$end] <> $qid\[$qb-$qe]\n";
        $qb > 0 ||  die "error: $tid\[$beg-$end] <> $qid\[$qb-$qe]\n";
        print $fho join("\t", $qid, $qb, $qe, $qsrd, $cid)."\n";
      }
      $cid ++;
    }
  }
  close $fhi;
  close $fho;
}
sub ucstat {
  my ($fi) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  my $h;
  my ($nc, $ns) = (0, 0);
  my ($lt, $li, $ln, $lc) = (0, 0, 0, 0);
  while(<$fhi>) {
    chomp;
    my ($cid, $id, $beg, $end, $str) = split "\t";
    my $len = $end - $beg + 1;
    $lt += $len;
    $h->{$id} ||= 1;
    if(defined($str) && $str ne "") {
      my @ps = split(" ", $str);
      $ln += $len;
      $lc += (@ps + 1) * $len;
      $lt += @ps * $len;
      for (@ps) {
        my @pps = split /\|/;
        $h->{$pps[0]} ||= 1;
      }
    } else {
      $li += $len;
    }
    $nc = max($nc, $cid);
  }
  $ns = scalar(keys(%$h));
  my $redcy = $lc / $lt;
  printf "%6d seqs, %6d clusters\n", $ns, $nc;
  printf "total [%10d]: sgnt [%10d] + clsr    [%10d]\n", $lt, $li, $lc;
  printf "dedup [%10d]: sgnt [%10d] + clsr-nr [%10d]\n", $li+$ln, $li, $ln;
  printf "redcy: %.05f\n", $lc / $lt;
  $lt == $li + $lc || die "len error\n";
  return $redcy;
}
sub clu2bed {
  my ($fi, $fo, $minlen) = @_;
  my ($hl, $hc) = read_cluster($fi);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  for my $cid (sort(keys(%$hc))) {
    my $ref = $hc->{$cid};
    my ($sid, $idx) = @{$ref->[0]};
    my ($beg, $end, $srd, $cid2) = @{$hl->{$sid}->[$idx]};
    $cid == $cid2 || die "error clu: $cid : $sid:$beg-$end\n";
    my $len = $end - $beg + 1;
    next if $len < $minlen;
    print $fho join("\t", $sid, $beg - 1, $end)."\n";
  }
  close $fho;
  runCmd("mv $fo tmp");
  runCmd("sortBed -sizeD -i tmp > $fo");
  runCmd("rm tmp");
}
sub clu2aln {
  my ($fi, $fs, $fo) = @_;
  my $db = Bio::DB::Fasta->new($fs);
  my ($hl, $hc) = read_cluster($fi);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  for my $cid (sort(keys(%$hc))) {
    my $ref = $hc->{$cid};
    next if @$ref <= 1;
    my (@strs, @seqs);
    for (@$ref) {
      my ($sid, $idx) = @$_;
      my ($beg, $end, $srd, $cid2) = @{$hl->{$sid}->[$idx]};
      $cid == $cid2 || die "error clu: $cid : $sid:$beg-$end\n";
      
      my $seq = $db->seq($sid, $beg, $end);
      $seq = Bio::Seq->new(-seq=>$seq)->revcom->seq if $srd eq "-";
      push @strs, "$sid-$beg-$end";
      push @seqs, $seq;
    }
    print $fho ">$cid ".join(" ", @strs)."\n";
    print $fho join("\n", @seqs)."\n";
  }
  close $fho;
}
sub get_fas_len {
  my ($fi) = @_;
  my $seqh = Bio::SeqIO->new(-file => "<$fi", -format=>'fasta');
  my $len = 0;
  while(my $seq = $seqh->next_seq()) {
    $len += $seq->length;
  }
  return $len;
}

