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

my $dir = "$ENV{'misc3'}/comp.ortho";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $tgt = "HM101";
my @qrys = qw/
  HM058 HM125 HM056 HM129 HM060
  HM095 HM185 HM034 HM004 HM050 
  HM023 HM010 HM022 HM324 HM340
/;
my @orgs = ($tgt, @qrys);

### run comp.syn.ortho.pl generate HM*_HM101/51_ortho/01.ortho.tbl
### run comp.ortho.ins.R create 01.tbl, 05.ortho.tbl, 06.no.ortho.tbl

#get_seq_by_org("05.ortho.tbl", "05.fas");
#get_seq_by_org("06.no.ortho.tbl", "06.fas");
#parallel_blat_8("06.fas", "05.fas", "08_blat", "08.blat.psl");
#runCmd("psl2gal.pl -i 08.blat.psl | gal.filter.ortho.pl | gal2tbl.pl -o 09.tbl");
#runCmd("sort -k1,1 09.tbl -o 09.tbl");
### run comp.ortho.R create 11.*.tbl 15.no.ortho.tbl


#get_seq_by_org("15.no.ortho.tbl", "15.fas");
#parallel_blat_8("15.fas", "15.fas", "16_blat", "16.blat.psl");
#runCmd("psl2gal.pl -i 16.blat.psl | gal.filter.ortho.pl | gal2mcl.pl -o 16.tbl");
#runCmd("\$soft/mcl/bin/mcl 16.tbl -te 4 -I 5.0 --abc -o 16.mcl");
#parse_mcl("16.mcl", "17.group.tbl", \@orgs);
### run comp.ortho.R create 21.*.tbl

ortho_group_cat("21.gid.tbl", "22.cat.tbl");
#ortho_aln_prep("21.gid.tbl", "25.aln.cmd", "25_seq", "25_aln");
#runCmd("parallel -j 16 --no-notice < 25.aln.cmd");

### run comp.ortho.R generate 28.dist.tbl
##get_seq_by_org("29.singleton.tbl", "29.fas");

#pull_orph_loc("31.ortho.tbl", "36.loc", "36.loc.cmd", \@qrys, $tgt);
#runCmd("parallel -j 16 --no-notice < 36.loc.cmd");
#pull_orph_loc_2("31.ortho.tbl", "36.loc", "36.loc.tbl");

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
      my $nid = defined($idx) ? "$org-$gid-$idx" : "$org-$gid";
      $seqHO->write_seq(Bio::Seq->new(-id => $nid, -seq => $seq));
    }
  }
  $seqHO->close();
}
sub parallel_blat_8 {
  my ($fi, $db, $do, $fo) = @_;
  runCmd("qsub.blat.pl -i $fi -o $do");
  runCmd("seq 0 7 | xargs -i printf \"%02d\\n\" {} | parallel --no-notice -j 8 blat -prot $db $do/part.{}.fas $do/part.{}.tbl");
  runCmd("seq 8 15 | xargs -i printf \"%02d\\n\" {} | parallel --no-notice -j 8 blat -prot $db $do/part.{}.fas $do/part.{}.tbl");
  runCmd("cat $do/part.*.tbl > $fo");
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
  my ($fi, $fo, $orgs) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", @$orgs)."\n";
  while(<$fhi>) {
    chomp;
    my @ps = split "\t";
    my $h = { map {$_ => []} @$orgs };
    my $hn = { map {$_ => 0} @$orgs };
    for (@ps) {
      my ($org, $id) = split /\-/;
      push @{$h->{$org}}, $id;
      $hn->{$org} ++;
    }
    my $n = max(values(%$hn));
    for my $i (1..$n) {
      my @ort;
      for my $org (@$orgs) {
        if($hn->{$org} >= $i) {
          push @ort, $h->{$org}->[$i-1];
        } else {
          push @ort, '';
        }
      }
      print $fho join("\t", @ort)."\n";
    }
  }
  close $fhi;
  close $fho;
}
sub ortho_group_refine {
  my ($fi, $fr, $fo, $orgs) = @_;
  my %hi = map {$_ => $orgs->[$_]} (0..@$orgs-1);
  my $ti = readTable(-in => $fi, -header => 1);
  $ti->delCols(['cat', 'n_org']);

  my $tr = readTable(-in => $fr, -header => 1);
  my $hr;
  for my $i (0..$tr->lastRow) {
    my @ps = $tr->row($i);
    if($ps[0] ne '') {
      exists $hr->{$ps[0]} && die "$ps[0] appeared >1 times in $fr\n";
      $hr->{$ps[0]} = \@ps;
    }
  }

  my $nc = 0;
  for my $i (0..$ti->lastRow) {
    my @ps = $ti->row($i);
    my $id = $ps[0];
    exists $hr->{$id} || next;
    my @psr = @{$hr->{$id}};
    my @psn = ('') x @$orgs;
    for my $j (1..@$orgs-1) {
      if($ps[$j] eq 'NA') {
        $ps[$j] = $psr[$j] if $psr[$j] ne '';
      } else {
        if($psr[$j] ne '' && $psr[$j] ne $ps[$j]) {
          $psn[$j] = $psr[$j];
        }
      }
    }
    my $nd = scalar( grep {$_ ne ''} @psn );
    if($nd > 0) {
      $nc ++;
      $ti->addRow(\@psn);
    }
  }
  print "$nc\n";
  
  for my $i (0..$tr->lastRow) {
    my @ps = $tr->row($i);
    $ti->addRow(\@ps) if $ps[0] eq '';
  }

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $ti->tsv(1);
  close $fho;
}
sub ortho_group_cat {
  my ($fi, $fo) = @_;
  my $ti = readTable(-in => $fi, -header => 1);
  my @cnames = $ti->header;
  my @orgs = @cnames[1..$#cnames];
  
  my $h;
  for my $org (@orgs) {
    my $fg = "$ENV{'genome'}/$org/51.gtb";
    my $tg = readTable(-in => $fg, -header => 1);
    my %hg = map {$tg->elm($_, 'id') => 
      [$tg->elm($_, 'cat2'), $tg->elm($_, 'cat3')]} (0..$tg->lastRow);
    $h->{$org} = \%hg;
  }
  
  my (@cats2, @cats3);
  for my $i (0..$ti->lastRow) {
    my ($idx, @gids) = $ti->row($i);
    my $hc;
    for my $j (0..$#gids) {
      my $gid = $gids[$j];
      my $org = $orgs[$j];
      if($gid ne '' && $gid ne '-') {
        my ($cat2, $cat3);
        if(!exists $h->{$org}->{$gid}) {
          ($cat2, $cat3) = ("Unknown", "");
        } else {
          ($cat2, $cat3) = @{$h->{$org}->{$gid}};
        }
        die "$org $gid\n" if !$cat2;
        $hc->{$cat2} ||= [0, ''];
        $hc->{$cat2}->[0] ++;
        $hc->{$cat2}->[1] = $cat3;
      }
    }
    my @cs2 = sort {$hc->{$b}->[0] <=> $hc->{$a}->[0]} keys(%$hc);
    my $cat2 = $cs2[0];
    my $cat3 = $hc->{$cat2}->[1];
    push @cats2, $cat2;
    push @cats3, $cat3;
  }
  $ti->addCol(\@cats2, 'cat2');
  $ti->addCol(\@cats3, 'cat3');
  $ti->delCols(\@orgs);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $ti->tsv(1);
  close $fho;
}
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
sub pull_orph_loc {
  my ($fi, $do, $fo, $qrys, $tgt) = @_;
  my $ti = readTable(-in => $fi, -header => 1);
  
  my ($hg, $hc);
  for my $qry (@$qrys) {
    my $fg = "$ENV{'genome'}/$qry/51.tbl";
    my $tg = readTable(-in => $fg, -header => 0);
    for my $i (0..$tg->lastRow) {
      my ($chr, $beg, $end, $srd, $id, $type) = $tg->row($i);
      $type eq "mrna" || next;
      $hg->{$qry}->{$id} = [$chr, $beg, $end];
    }
    my $fgax = "$ENV{'misc3'}/$qry\_$tgt/23_blat/41.9/gax.gz";
    $hc->{$qry} = Tabix->new(-data => $fgax);
  }

  my $t0 = [gettimeofday];
  my @orgs = $ti->header;
  my @idxs_orph = indexes {$_ eq ""} $ti->col($tgt);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  for my $i (1..@idxs_orph) {
    my $idx = $idxs_orph[$i-1];
    my @idsa = $ti->row($idx);
    $idsa[0] eq "" || die "$idx has $tgt id\n";

    my @beds;
    my @js = indexes {$_ ne ""} @idsa;
    for my $j (@js) {
      my ($org, $id) = ($orgs[$j], $idsa[$j]);
      exists $hg->{$org}->{$id} || die "$org: $id not found\n";
      my ($c, $b, $e) = @{$hg->{$org}->{$id}};
      my $ary = read_gax($hc->{$org}, $c, $b, $e);
      @$ary == 0 && next;
      for (@$ary) {
        my ($cid, $tid, $tb, $te, $tsrd, $qid, $qb, $qe) = @$_;
        push @beds, [$qid, $qb-1, $qe, 1];
      }
    }
    @beds > 0 || next;
    my $ft = "$do/$idx.bed";
    open(my $fht, ">$ft") or die "cannot write $ft\n";
    for (@beds) { print $fht join("\t", @$_)."\n" }
    close $fht;

    print $fho "sortBed -i $ft | mergeBed -i stdin -c 4 -o sum | sort -k4,4nr | head -n 1 | cut -f1-3 > $do/$idx.out\n";
    my $t1 = [gettimeofday];
    printf "***%.01f min $i\n", tv_interval($t0, $t1) / 60 if $i % 500 == 0;
  }
  close $fho;
}
sub pull_orph_loc_2 {
  my ($fi, $di, $fo) = @_;
  my $ti = readTable(-in => $fi, -header => 1);
  
  my @orgs = $ti->header;
  my @idxs_orph = indexes {$_ eq ""} $ti->col($tgt);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  for my $i (1..@idxs_orph) {
    my $idx = $idxs_orph[$i-1];

    my $ft = "$di/$idx.out";
    -s $ft || next;
    open(my $fht, "<$ft") or die "cannot read $ft\n";
    my $line = <$fht>;
    close $fht;
    print $fho join("\t", $idx, $line);
  }
  close $fho;
}
__END__

