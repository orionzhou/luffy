#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  pan4probe.pl 

=head1 SYNOPSIS
  
  pan4probe.pl [-help]

  Options:
    -h (--help)   brief help message
    -s (--stat)   print statistics

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Tabix;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Data::Dumper;
use Common;
use Location;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $help_flag;
my $stat_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "stat|s"  => \$stat_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $dir = "$ENV{'misc3'}/pan4probe";
chdir $dir || die "cannot chdir to $dir\n";

##### design back-bone probes
#runCmd("probedesign -i 03.uniq.bed -d 500 -o 06.prb");
#agilentbed2prb('90_agilent_pan3/Tiling_bed', '06.prb');

##### design gene-enriched probes
#prepare_genic("$ENV{'genome'}/HM101/55_noTE.gtb", "11.genic.tbl");
#prepare_novel("$ENV{'misc3'}/pan4seq/33.agp.tbl", "12.novel.tbl");
#runCmd("cat 11.genic.tbl 12.novel.tbl > 15.tbl");
#tbl2bed_genic("15.tbl", "15.bed");
#runCmd("sortBed -i 15.bed | mergeBed -i stdin > 15.merged.bed");

### run R script - generate partitioned Beds
#design_probes("15.tbl", "21.cds.nov", "22.prb", 4);
#runCmd("ln -sf 21.cds.nov/rd4/26.prb 22.prb");

##### merge two probe groups
#merge_backbone_genic("06.prb", "22.prb");

#stat_probe("15.tbl", "36.prb.gz", "41.tbl");

#prb2tdt('36.prb', '61.tdt', "$ENV{'genome'}/pan4/11_genome.fas");

sub agilentbed2prb {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot write $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  while(<$fhi>) {
    chomp;
    next unless /^chr/i;
    my ($chr, $beg, $end) = split "\t";
    print $fho join("\t", $chr, $beg + 1)."\n";
  }
  close $fhi;
  close $fho;
}
sub prepare_genic {
  my ($fi, $fo) = @_;
  my $ti = readTable(-in => $fi, -header => 1);

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/chr beg end id type loc/)."\n";
  for my $i (0..$ti->lastRow) {
    my ($id, $par, $chr, $beg, $end, $srd, 
      $locES, $locIS, $locCS, $loc5S, $loc3S, $phase, 
      $src, $conf, $cat1, $cat2, $cat3, $note) = $ti->row($i);
    $cat3 = "CRP" if $cat3 =~ /crp/i;
    $cat3 = "NBS" if $cat3 =~ /nbs/i;
    my $rloc = locStr2Ary($locCS);
    my $aloc = $srd eq "+" ? 
      [ map {[$beg + $_->[0] - 1, $beg + $_->[1] - 1]} @$rloc ] :
      [ map {[$end - $_->[1] + 1, $end - $_->[0] + 1]} reverse @$rloc ];
    my $loc = locAry2Str($aloc);
    print $fho join("\t", $chr, $beg, $end, $id, $cat3, $loc)."\n";
  }
  close $fho;
}
sub prepare_novel {
  my ($fi, $fo) = @_;

  my $h;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  while(<$fhi>) {
    chomp;
    next if /^cid/;
    my @ps = split "\t";
    my ($chr, $beg, $end) = @ps[11..13];
    my $loc = "$beg-$end";
    $h->{$chr} ||= 0;
    $h->{$chr} ++;
    my $id = sprintf "%s.%06d", $chr, $h->{$chr};
    print $fho join("\t", $chr, $beg, $end, $id, "pchr", $loc)."\n";
  }
  close $fhi;
  close $fho;
}
sub tbl2bed_genic {
  my ($fi, $fo) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  for my $i (0..$t->lastRow) {
    my ($chr, $beg, $end, $id, $type, $locS) = $t->row($i);
    my $loc = locStr2Ary($locS);
    for (@$loc) {
      my ($b, $e) = @$_;
      print $fho join("\t", $chr, $b - 1, $e)."\n";
    }
  }
  close $fho;
}
sub tbl2bed_failed {
  my ($fi, $fo) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  for my $i (0..$t->lastRow) {
    my ($chr, $beg, $end, $id, $type, $locS, $prbnum) = $t->row($i);
    $prbnum < 3 || next;
    my $loc = locStr2Ary($locS);
    for (@$loc) {
      my ($b, $e) = @$_;
      print $fho join("\t", $chr, $b - 1, $e)."\n";
    }
  }
  close $fho;
}
sub read_prb {
  my ($con, $id, $beg, $end) = @_;
  my $iter = $con->query($id, $beg - 1, $end);
  my @ary;
  $iter->get || return \@ary;
  while (my $line = $con->read($iter)) {
    my @ps = split("\t", $line);
    @ps > 0 || next;
    my ($chr, $pos) = @ps;
    push @ary, $pos;
  }
  return \@ary;
}
sub read_prb_by_loc {
  my ($con, $id, $loc) = @_;
  my @ary;
  for (@$loc) {
    my ($beg, $end) = @$_;
    my $sry = read_prb($con, $id, $beg, $end);
    push @ary, @$sry if @$sry > 0;
  }
  return \@ary;
}
sub stat_probe {
  my ($fi, $fp, $fo) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  my $con = Tabix->new(-data => $fp);

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/chr beg end id type loc prbnum prbpos/)."\n";

  my (@prbnum, @prbpos);
  for my $i (0..$t->lastRow) {
    my ($chr, $beg, $end, $id, $type, $locS) = $t->row($i);
    my $loc = locStr2Ary($locS);
    my $poss = read_prb_by_loc($con, $chr, $loc);
    print $fho join("\t", $chr, $beg, $end, $id, $type, $locS,
      scalar(@$poss), join(",", @$poss))."\n";
  }
  close $fho;
}
sub stat_probe_two {
  my ($fi, $fp1, $fp2, $fo) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  my $con1 = Tabix->new(-data => $fp1);
  my $con2 = Tabix->new(-data => $fp2);

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/chr beg end id type loc 
    prbnum1 prbpos1 prbnum2 prbpos2/)."\n";

  my (@prbnum1, @prbpos1, @prbnum2, @prbpos2);
  for my $i (0..$t->lastRow) {
    my ($chr, $beg, $end, $id, $type, $locS) = $t->row($i);
    my $loc = locStr2Ary($locS);
    my $prbpos1 = read_prb_by_loc($con1, $chr, $loc);
    my $prbpos2 = read_prb_by_loc($con2, $chr, $loc);
    print $fho join("\t", $chr, $beg, $end, $id, $type, $locS,
      scalar(@$prbpos1), join(",", @$prbpos1),
      scalar(@$prbpos2), join(",", @$prbpos2))."\n";
  }
  close $fho;
}
sub get_equal_idx {
  my ($n, $a) = @_;
  $n <= $a || die "cannot select $n out of $a\n";
  my $dist = $a / $n;
  return map { int($dist * $_) } (0..$n-1);
}
sub pick_more_probe {
  my ($fi, $fp, $n_prb) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  open(my $fhp, ">$fp") or die "cannot write $fp\n";
  for my $i (0..$t->lastRow) {
    my ($chr, $beg, $end, $id, $type, $locS, $prbnum1, $prbpos1, 
      $prbnum2, $prbpos2) = $t->row($i);
    my $n = min( $n_prb - min($n_prb, $prbnum1), $prbnum2 );
    if($n > 0) {
      my @poss1 = split(",", $prbpos1);
      my @poss2 = split(",", $prbpos2);
      my @idxs = get_equal_idx($n, $prbnum2);
      for (@idxs) {
        print $fhp join("\t", $chr, $poss2[$_])."\n";
      }
      my @poss = sort {$a <=> $b} (@poss1, @poss2[@idxs]);
    }
  }
  close $fhp;
}
sub fix_prbstat_rd1 {
  my ($fi) = @_;
  runCmd("mv $fi $fi.tmp");
  my $t = readTable(-in => "$fi.tmp", -header => 1);

  my $prbnum1 = [(0) x $t->nofRow];
  my $prbpos1 = [('') x $t->nofRow];
  my $prbnum2 = $t->delCol("prbnum");
  my $prbpos2 = $t->delCol("prbpos");
  $t->addCol($prbnum1, 'prbnum1');
  $t->addCol($prbpos1, 'prbpos1');
  $t->addCol($prbnum2, 'prbnum2');
  $t->addCol($prbpos2, 'prbpos2');
  open(my $fho, ">$fi") or die "cannot write to $fi\n";
  print $fho $t->tsv(1);
  close $fho;

  runCmd("rm $fi.tmp");
}
sub design_probes {
  my ($fi, $dir, $fo, $n_prb) = @_;
  chdir $dir || die "cannot chdir to $dir\n";
  
  my $rd = 6;
  for my $i (1..$rd) {
    print "########## working on rd $i\n";
#    $i > 1 || next;
    my $dw = "rd$i";
    -d $dw || make_path($dw);
    chdir $dw || die "cannot chdir to $dw\n";
    
    my $fb = "../rd$i.bed";
    my $dp = "../rd".($i - 1);
    if($i == 1) {
      runCmd("ln -sf $fb 10.bed");
    } else {  
      runCmd("subtractBed -a $fb -b $dp/28.prbexpand.bed > 10.bed");
    }
    runCmd("probedesign.pl -i 10.bed -d 150 -o 11.prb");
    runCmd("bgzip -c 11.prb > 11.prb.gz");
    runCmd("tabix -s 1 -b 2 -e 2 11.prb.gz");
    if($i == 1) {
      stat_probe("../../$fi", "11.prb.gz", "21.tbl");
      fix_prbstat_rd1("21.tbl");
      pick_more_probe("21.tbl", "22.prb", $n_prb);
      runCmd("cp 22.prb 26.prb");
    } else {
      stat_probe_two("$dp/31.tbl", "$dp/26.prb.gz", "11.prb.gz", "21.tbl");
      pick_more_probe("21.tbl", "22.prb", $n_prb);
      runCmd("cat $dp/26.prb 22.prb > 26.prb");
    }
    runCmd("sort -k1,1 -k2,2n -u 26.prb -o 26.prb");
    runCmd("bgzip -c 26.prb > 26.prb.gz");
    runCmd("tabix -s 1 -b 2 -e 2 26.prb.gz");
    runCmd("probeexpand.pl -i 26.prb -d 150 -o 28.prbexpand.bed");
    stat_probe("../../$fi", "26.prb.gz", "31.tbl");

    chdir "..";
  }
  chdir "..";
}
sub prb2bb {
  my ($fi, $fo, $fs) = @_;
  runCmd("awk 'BEGIN{OFS=\"\\t\"} {print \$1, \$2-1, \$2+59}' $fi > $fi.bed");
  runCmd("bedToBigBed -tab $fi.bed $fs $fo");
  runCmd("rm $fi.bed");
}
sub merge_backbone_genic {
  my ($fb, $fg) = @_;
  runCmd("probeexpand.pl -i $fb -d 0 -o 31.backbone.bed");
  runCmd("probeexpand.pl -i $fg -d 150 -o 32.genic.exp.bed");
  runCmd("subtractBed -A -a 31.backbone.bed -b 32.genic.exp.bed > 35.backbone.flt.bed");
  runCmd("cut -f1,3 35.backbone.flt.bed > 35.backbone.flt.prb");
  runCmd("cat 22.prb 35.backbone.flt.prb > 36.prb");
  
  runCmd("sort -k1,1 -k2,2n -u 36.prb -o 36.prb");
  runCmd("bgzip -c 36.prb > 36.prb.gz");
  runCmd("tabix -s 1 -b 2 -e 2 36.prb.gz");
  prb2bb("36.prb", "36.bb", "$ENV{'genome'}/pan4/15.sizes");
}
sub prb2tdt {
  my ($fi, $fo, $fs) = @_;
  my $db = Bio::DB::Fasta->new($fs);
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/ProbeID Sequence TargetID Accessions GeneSymbols Description ChromosomalLocation/)."\n";
  my $cnt = 1;
  while(<$fhi>) {
    chomp;
    my ($chr, $beg) = split "\t";
    my $end = $beg + 60 - 1;
    my $seq = $db->seq($chr, $beg, $end);
    my $id = sprintf "mtpan4_%06d", $cnt++;
    print $fho join("\t", $id, $seq, ('') x 4, "$chr:$beg-$end")."\n";
  }
  close $fhi;
  close $fho;
}

