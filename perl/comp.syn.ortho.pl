#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.ortho.syn.pl - find orthologs in pairwise comparison

=head1 SYNOPSIS
  
  comp.ortho.syn.pl [-help] [-qry qry-genome] [-tgt tgt-genome]

  Options:
    -h (--help)   brief help message
    -q (--qry)    qry genome
    -t (--tgt)    tgt genome

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
use Tabix;
use Bio::DB::Fasta;
use Data::Dumper;
use Common;
use Location;
use Gtb;
use Gal;
use List::Util qw/min max sum/;

my ($qry, $tgt) = qw/HM034 HM101/;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "qry|q=s"  => \$qry,
  "tgt|t=s"  => \$tgt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$qry || !$tgt;

my $dir = "$ENV{'misc3'}/$qry\_$tgt/51_ortho";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $qdir = "$ENV{'genome'}/$qry";
my $tdir = "$ENV{'genome'}/$tgt";
my $cdir = "$ENV{'misc3'}/$qry\_$tgt/23_blat";

my $gax = Tabix->new(-data => "$cdir/31.9/gax.gz");
my $snp = Tabix->new(-data => "$cdir/31.9/snp.gz");

my $fgq = "$qdir/51.gtb";
my $fgt = "$tdir/51.gtb";
my $qdb = Bio::DB::Fasta->new("$qdir/51.fas");
my $tdb = Bio::DB::Fasta->new("$tdir/51.fas");

my $qgene = Tabix->new(-data => "$qdir/51.tbl.gz");

syn_ortho($fgq, $fgt, $gax, $qgene, "01.ortho.tbl");
#runCmd("blat -prot $tdir/51.fas $qdir/51.fas 11.psl");
#runCmd("psl2gal.pl -i 11.psl | gal.filter.ortho.pl -o 12.gal");
#syn_ortho_refine("01.ortho.tbl", "12.gal", "21.ortho.tbl");
#ortho_combine("21.ortho.tbl", $fgq, $fgt,  "12.gal", "31.ortho.tbl");

sub syn_ortho {
  my ($fgq, $fgt, $gax, $qgene, $fo) = @_;
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/idx tid tlen qid qlen slen cid lev qstr/)."\n";

  my $tq = readTable(-in => $fgq, -header => 1);
  my $hq;
  for my $i (0..$tq->lastRow) {
    my ($qid, $qpar, $qchr, $qb, $qe, $qsrd, 
      $elocS, $ilocS, $clocS, $flocS, $tlocS, $phase) = $tq->row($i);
    my $qloc = locStr2Ary($clocS);
    $hq->{$qid} = ['', locAryLen($qloc)];
  }

  my $tt = readTable(-in => $fgt, -header => 1);
  for my $i (1..$tt->nofRow) {
    my ($tid, $tpar, $tchr, $tb, $te, $tsrd, 
      $elocS, $ilocS, $clocS, $flocS, $tlocS, $phase, 
      $src, $conf, $cat1, $cat2, $cat3, $note) = $tt->row($i - 1);
    my $tloc = locStr2Ary($clocS); 
    $tloc = $tsrd eq "+" ? [ map {[$tb+$_->[0]-1, $tb+$_->[1]-1]} @$tloc ]
      : [ map {[$te-$_->[1]+1, $te-$_->[0]+1]} @$tloc ];
    my $tlen = locAryLen($tloc);
    
    my $h;
    for (@$tloc) {
      my ($tbeg, $tend) = @$_;
      my $ary = read_gax($gax, $tchr, $tbeg, $tend);
      for (@$ary) {
        my ($cid, $tid1, $tb1, $te1, $tsrd1, $qid, $qb, $qe, $qsrd, $lev) = @$_;
        my $ary2 = read_cds($qgene, $qid, $qb, $qe);
        for (@$ary2) {
          my ($chr, $beg, $end, $srd, $gene) = @$_;
          $h->{$gene} ||= [0, $cid, $lev];
          $h->{$gene}->[0] += $end - $beg + 1;
        }
      }
    }
    if(!defined($h)) {
      print $fho join("\t", $i, $tid, $tlen, ('') x 6)."\n";
    } else {
      my @qids = sort {$h->{$b}->[0] <=> $h->{$a}->[0]} keys(%$h);
      my $qid = $qids[0];
      my ($slen, $cid, $lev) = @{$h->{$qid}};
      my $qlen = $hq->{$qid}->[1];
      my $str = join(",", map {$_."/".$hq->{$_}->[1]."/".$h->{$_}->[0]} keys(%$h));
      print $fho join("\t", $i, $tid, $tlen, $qid, $qlen, $slen, $cid, $lev, $str)."\n";
    }
  }
  close $fho;
}
sub syn_ortho_raw {
  my ($fgq, $fgt, $gax, $qgene, $fo) = @_;

  my $tq = readTable(-in => $fgq, -header => 1);
  my $hq;
  for my $i (0..$tq->lastRow) {
    my ($qid, $qpar, $qchr, $qb, $qe, $qsrd, 
      $elocS, $ilocS, $clocS, $flocS, $tlocS, $phase) = $tq->row($i);
    my $qloc = locStr2Ary($clocS);
    $hq->{$qid} = ['', locAryLen($qloc)];
  }

  my $tt = readTable(-in => $fgt, -header => 1);
  my $ht;
  for my $i (0..$tt->lastRow) {
    my ($tid, $tpar, $tchr, $tb, $te, $tsrd, 
      $elocS, $ilocS, $clocS, $flocS, $tlocS, $phase, 
      $src, $conf, $cat1, $cat2, $cat3, $note) = $tt->row($i);
    my $tloc = locStr2Ary($clocS); 
    $tloc = $tsrd eq "+" ? [ map {[$tb+$_->[0]-1, $tb+$_->[1]-1]} @$tloc ]
      : [ map {[$te-$_->[1]+1, $te-$_->[0]+1]} @$tloc ];
    my $tlen = locAryLen($tloc);
    
    my $h;
    for (@$tloc) {
      my ($tbeg, $tend) = @$_;
      my $ary = read_gax($gax, $tchr, $tbeg, $tend);
      for (@$ary) {
        my ($cid, $tid1, $tb1, $te1, $tsrd1, $qid, $qb, $qe, $qsrd) = @$_;
        my $ary2 = read_cds($qgene, $qid, $qb, $qe);
        for (@$ary2) {
          my ($chr, $beg, $end, $srd, $gene) = @$_;
          $h->{$gene} ||= 0;
          $h->{$gene} += $end - $beg + 1;
        }
      }
    }
    if(!defined($h)) {
      $ht->{$tid} = '';
      next;
    }
    my @qids = sort {$h->{$b} <=> $h->{$a}} keys(%$h);
    my $qid = $qids[0];
    my $qlen = $hq->{$qid}->[1];
    my $slen = $h->{$qid};
    if($slen / $tlen >= 0.2 && $slen / $qlen >= 0.2) {
#      my $tseq = $tdb->seq($tid);
#      my $qseq = $qdb->seq($qid);
#      $tseq =~ s/X$//i;
#      $qseq =~ s/X$//i;
      if($hq->{$qid}->[0] ne "") {
        my $ptid = $hq->{$qid}->[0];
        if($ht->{$ptid}->[3] >= $slen) {
          $ht->{$tid} = '';
          next;
        }
        $ht->{$tid} = [$tlen, $qid, $qlen, $slen];
        $hq->{$qid}->[0] = $tid;
        $ht->{$ptid} = '';
      } else {
        $ht->{$tid} = [$tlen, $qid, $qlen, $slen];
        $hq->{$qid}->[0] = $tid;
      }
    } else {
      $ht->{$tid} = '';
    }
  }

  my ($nt, $nq) = ($tt->nofRow, $tq->nofRow);
  my ($no, $nst, $nsq) = (0) x 3;

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/tid tlen qid qlen slen/)."\n";
  for my $tid (sort(keys(%$ht))) {
    if($ht->{$tid} eq "") {
      $nst ++;
    } else {
      my $ps = $ht->{$tid};
      print $fho join("\t", $tid, @$ps)."\n";
      $no ++;
    }
  }
  close $fho;
 
  for my $qid (sort(keys(%$hq))) {
    $hq->{$qid}->[0] eq "" || next;
    $nsq ++;
  }

  $nt == $no + $nst || die "error: nt[$nt] <> no[$no] + nst[$nst]\n";
  $nq == $no + $nsq || die "error: nq[$nq] <> no[$no] + nst[$nsq]\n";
  printf "tgt|%5d: syn-ortho|%5d + singleton|%5d\n", $nt, $no, $nst;
  printf "qry|%5d: syn-ortho|%5d + singleton|%5d\n", $nq, $no, $nsq;
}
sub needle_ident {
  my ($seq1, $seq2) = @_;
  runCmd("echo -e \">seq1\\n$seq1\" > seq1.fas", 0); 
  runCmd("echo -e \">seq2\\n$seq2\" > seq2.fas", 0); 
  runCmd("needle seq1.fas seq2.fas -gapopen 10 -gapextend 0.5 \\
    -aformat pair -outfile aln.needle", 0);
  open(my $fha, "<aln.needle") or die "cannot read aln.needle\n";
  my ($len, $idt, $sim, $gap, $sco);
  while(<$fha>) {
    chomp;
    $len = $1 if /^\# Length:\s*(\d+)/;
    $idt = $1 if /^\# Identity:\s*(\d+)/;
    $sim = $1 if /^\# Similarity:\s*(\d+)/;
    $gap = $1 if /^\# Gaps:\s*(\d+)/;
    $sco = $1 if /^\# Score:\s*([\d\.]+)/;
  }
  close $fha;
  runCmd("rm seq1.fas seq2.fas aln.needle", 0);
  my $ident = sprintf "%.03f", $idt / $len;
  return $ident;
}
sub read_cds {
  my ($con, $chr, $beg, $end) = @_;
  my $iter = $con->query($chr, $beg - 1, $end);
  my @ary;
  return \@ary if ! $iter->get();
  while (my $line = $con->read($iter)) {
    my @ps = split("\t", $line);
    my ($chr, $beg1, $end1, $srd, $id, $type, $cat) = @ps;
    $type eq "cds" || next;
    push @ary, [$chr, max($beg, $beg1), min($end, $end1), $srd, $id];
  }
  return \@ary;
}
sub syn_ortho_refine {
  my ($fi, $fb, $fo) = @_;
  my $hb;
  open(my $fhb, "<$fb") or die "cannot read $fb\n";
  while(<$fhb>) {
    chomp;
    next if /^id/;
    my @ps = split "\t";
    my ($tid, $qid, $score) = @ps[1,6,18];
    $hb->{$qid}->{$tid} = $score;
  }
  close $fhb;

  my @idxs;
  my $ti = readTable(-in => $fi, -header => 1);
  for my $i (0..$ti->lastRow) {
    my ($tid, $tlen, $qid, $qlen ,$slen) = $ti->row($i);
    push @idxs, $i if !exists $hb->{$qid}->{$tid};
  }
  
  $ti->delRows(\@idxs);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $ti->tsv(1);
  close $fho;
  my ($nr, $np) = (scalar(@idxs), $ti->nofRow);
  my $n = $nr + $np;
  printf "%5d | %5d passed, %5d removed\n", $np, $n, $nr; 
}
sub read_ids {
  my ($fi) = @_;
  my @ary;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  while(<$fhi>) {
    chomp;
    /(^\#)|(^\s*$)/ && next;
    push @ary, $_;
  }
  close $fhi;
  return \@ary;
}
sub ortho_combine {
  my ($fi, $fgq, $fgt, $fb, $fo) = @_;

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/tid qid/)."\n";

  my $ti = readTable(-in => $fi, -header => 1);;
  my ($hoq, $hot);
  for my $i (0..$ti->lastRow) {
    my ($tid, $tlen, $qid, $qlen ,$slen) = $ti->row($i);
    $hoq->{$qid} = $tid;
    $hot->{$tid} = [$tlen ,$qid, $qlen, $slen];
    print $fho join("\t", $tid, $qid)."\n";
  }

  my $hb;
  open(my $fhb, "<$fb") or die "cannot read $fb\n";
  while(<$fhb>) {
    chomp;
    next if /^id/;
    my @ps = split "\t";
    my ($tid, $qid, $score) = @ps[1,6,18];
    next if exists $hoq->{$qid} || $hot->{$tid};
    $hb->{$qid}->{$tid} = $score;
  }
  close $fhb;

  my $n = 0;
  my $hq;
  for my $qid (sort(keys(%$hb))) {
    my $hs = $hb->{$qid};
    my @tids = grep {!exists $hq->{$_}} keys(%$hs);
    next if @tids == 0;
    @tids = sort {$hs->{$b} <=> $hs->{$a}} @tids;
    my $tid = $tids[0];
    $hq->{$qid} = $tid;
    $n ++;
    print $fho join("\t", $tid, $qid)."\n";
  }
  close $fho;
  printf "syn-ortho|%5d + RBH-ortho|%5d = %5d\n", $ti->nofRow, $n, $n + $ti->nofRow;
}

__END__

