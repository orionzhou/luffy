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

my $dir = "$ENV{'misc3'}/comp.ortho.ins";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $tgt = "HM101";
my @qrys = qw/
  HM058 HM125 HM056 HM129 HM060
  HM095 HM185 HM034 HM004 HM050 
  HM023 HM010 HM022 HM324 HM340
/;
my @orgs = ($tgt, @qrys);

#make_ins_fas("01.tbl", "03_seq", \@qrys);
grouping_seq("01.tbl", "03_seq", "05.group.tbl", \@qrys);

sub make_ins_fas {
  my ($fi, $do, $qrys) = @_;
  -d $do || make_path($do);
  my $hdb;
  for my $org (@$qrys) {
    my $db = Bio::DB::Fasta->new("$ENV{'genome'}/$org/51.fas");
    $hdb->{$org} = $db;
  }
  my $ti = readTable(-in => $fi, -header => 1);
  for my $i (0..$ti->lastRow) {
    my ($idx, $chr, $pos, $note) = $ti->row($i);
    my @ary = split(" ", $note);
    @ary > 1 || next;

    my $id = sprintf("$do/%06d", $idx);
    open(my $fhs, ">$id.fas") or die "cannot write to $id.fas\n";
    for (@ary) {
      my ($org, $gid, $len) = split "-";
      my $seq = $hdb->{$org}->seq($gid);
      print $fhs join("\n", ">$org-$gid", $seq)."\n";
    }
    close $fhs;
#    runCmd("muscle -quiet -in $do/$kid.fas -out $do/$kid.aln.fas");
  }
}
sub grouping_seq {
  my ($fi, $di, $fo, $orgs) = @_;
  my $ti = readTable(-in => $fi, -header => 1);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/idx chr pos/, @$orgs)."\n";
  for my $i (0..$ti->lastRow) {
    my ($idx, $chr, $pos, $note) = $ti->row($i);
    my @ary = split(" ", $note);
    @ary > 1 || next;
    
#    next unless $idx == 14;
    my $fs = sprintf("$di/%06d.fas", $idx);
    runCmd("blat -prot $fs $fs x.psl");
    runCmd("psl2gal.pl -i x.psl | gal2mcl.pl -o x.txt");
    runCmd("\$soft/mcl/bin/mcl x.txt -te 4 -I 5.0 --abc -o x.mcl");
    my $groups = read_mcl("x.mcl");
    if(@$groups == 0) {
      for (@ary) {
        my ($org, $gid, $len) = split "-";
        push @$groups, ["$org-$gid"];
      }
    }
    my $ary = partition_record($groups, $orgs);
    for (@$ary) {
      my @orts = @$_;
      print $fho join("\t", $idx, $chr, $pos, @orts)."\n";
    }
  } 
  runCmd("rm x.psl x.txt x.mcl");
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
  my ($fi, $fo, $orgs) = @_;
  my $ti = readTable(-in => $fi, -header => 1);
  my $h;
  for my $org (@$orgs) {
    my $fg = "$ENV{'genome'}/$org/51.gtb";
    my $tg = readTable(-in => $fg, -header => 1);
    my %hg = map {$tg->elm($_, 'id') => 
      [$tg->elm($_, 'cat2'), $tg->elm($_, 'cat3')]} (0..$tg->lastRow);
    $h->{$org} = \%hg;
  }

  my (@cats2, @cats3);
  for my $i (0..$ti->lastRow) {
    my @ps = $ti->row($i);
    my $hc;
    for my $j (0..@$orgs-1) {
      if($ps[$j] ne '' && $ps[$j] ne 'NA') {
        my $org = $orgs->[$j];
        my $cat2 = $h->{$org}->{$ps[$j]}->[0];
        my $cat3 = $h->{$org}->{$ps[$j]}->[1];
        die "$org $ps[$j]\n" if !$cat2;
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
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $ti->tsv(1);
  close $fho;
}
sub ortho_aln_prep {
  my ($fi, $fo) = @_;
  -d "42.seq" || make_path("42.seq");
  -d "42.aln" || make_path("42.aln");
  my $ti = readTable(-in => $fi, -header => 1);
  my @orgs = $ti->header;

  my $h;
  for my $org (@orgs) {
    my $fs = "$ENV{'genome'}/$org/51.fas";
    $h->{$org} = Bio::DB::Fasta->new($fs);
  }

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  for my $i (0..$ti->lastRow) {
    my @ps = $ti->row($i);
    my (@ids, @seqs);
    for my $j (0..$#ps) {
      $ps[$j] eq '' && next;
      my $org = $orgs[$j];
      push @ids, "$org|$ps[$j]";
      push @seqs, $h->{$org}->seq($ps[$j]);
    }
    if(@ids > 1) {
      my $fn = sprintf("%06d", $i);
      open(my $fhs, ">42.seq/$fn.fas") or die "cannot write $fn\n";
      print $fhs join("\n", map {">".$ids[$_]."\n".$seqs[$_]} 0..$#ids);
      close $fhs;

      print $fho "clustalo -i 42.seq/$fn.fas -o 42.aln/$fn.fas --outfmt=fasta --force --full --full-iter\n";
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

