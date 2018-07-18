#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  pan4seq.pl - construct Pan4 genome fasta

=head1 SYNOPSIS
  
  pan4seq.pl [-help]

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
use Bio::SeqIO;
use Bio::DB::Fasta;
use Tabix;
use Data::Dumper;
use Common;
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

my $dir = "/home/youngn/zhoup/Data/misc3/pan4seq";
chdir $dir || die "cannot chdir to $dir\n";

my @orgs = qw/HM340 HM034 HM056/;

if($stat_flag) {
  cluster_stat();
  exit;
}

#merge_seq(\@orgs, '01.fas');
#runCmd("usearch.pl -i 01.fas -o 11");
#cluster_coord("11.clu", "13.clu");
#cluster_group("13.clu", "21.cluster.tbl", \@orgs);
#cluster_expand("21.cluster.tbl", "21.coord.tbl");
#cluster_dedup("21.cluster.tbl", "25.cluster.tbl");
#cluster_expand("25.cluster.tbl", "25.coord.tbl");
#cluster_filter("25.cluster.tbl", "27.cluster.tbl", 50);
#cluster_expand("27.cluster.tbl", "27.coord.tbl");

#cluster_filter("25.cluster.tbl", "28.cluster.tbl", 60);
#cluster_seq_prepare("28.cluster.tbl", "31.tbl", \@orgs);
#seq_coord_anchor("31.tbl", "32.anchored.tbl", \@orgs);
#cluster_seq("32.anchored.tbl", "33.agp.tbl", "33.pan4.fas");

sub merge_seq {
  my ($orgs, $fo) = @_;
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  my $seqHO = Bio::SeqIO->new(-fh=>$fho, -format=>'fasta');
  my ($len, $len2) = (0, 0);
  for my $org (@$orgs) {
    my $fd = "/home/youngn/zhoup/Data/genome/$org/11_genome.fas";
    my $db = Bio::DB::Fasta->new($fd);
    my $fi = "../$org\_HM101/41_novseq/41.bed";
    open(my $fhi, "<$fi") or die "cannot read $fi\n";
    while(<$fhi>) {
      chomp;
      my ($id, $beg, $end) = split "\t";
      $beg += 1;
      $len += $end - $beg + 1;
      my $nid = "$org-$id-$beg-$end";
      my $seq = $db->seq($id, $beg, $end);
      $seqHO->write_seq( Bio::Seq->new(-id=>$nid, -seq=>$seq) );
      if($end > $db->length($id)) {
        printf "$org $id end[$end] > len[%d]\n", $db->length($id);
      }
      $len2 += length($seq);
    }
    close $fhi;
  }
  close $fho;
  printf "bedlen: %8d\n", $len;
  printf "seqlen: %8d\n", $len2;
}
sub cluster_coord {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  while(<$fhi>) {
    chomp;
    my ($idstr, $rb, $re, $srd, $cid) = split "\t";
    my ($org, $id, $beg, $end) = split("-", $idstr);
    my ($b, $e) = ($beg + $rb - 1, $beg + $re - 1);
    print $fho join("\t", $org, $id, $b, $e, $srd, $cid)."\n";
  }
  close $fhi;
  close $fho;
}
sub cluster_group {
  my ($fi, $fo, $orgs_order) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  my $hc;
  while(<$fhi>) {
    chomp;
    my ($org, $id, $beg, $end, $srd, $cid) = split "\t";
    $hc->{$cid} ||= [];
    push @{$hc->{$cid}}, [$org, $id, $beg, $end, $srd];
  }
  close $fhi;

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/cid len cnt orgs strs/)."\n";
  for my $cid (sort {$a <=> $b} keys(%$hc)) {
    my (@strs, $ho);
    my $len;
    for (@{$hc->{$cid}}) {
      my ($org, $id, $beg, $end, $srd) = @$_;
      push @strs, "$org:$id:$beg:$end:$srd";

      $len ||= $end - $beg + 1;
      $len == $end - $beg + 1 || die "len conflict: cid $cid\n";
      $ho->{$org} = 1;
    }
    my $cnt = @{$hc->{$cid}};
    my @orgs = grep {exists $ho->{$_}} @$orgs_order;
    my $str_org = join(" ", @orgs);
    my $str_loc = join(" ", @strs);

    print $fho join("\t", $cid, $len, $cnt, $str_org, $str_loc)."\n";
  }
  close $fho;
}
sub cluster_expand {
  my ($fi, $fo) = @_;
  my $ti = readTable(-in => $fi, -header => 1);

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/cid org orgs id beg end srd/)."\n";
  for my $i (0..$ti->lastRow) {
    my ($cid, $len, $cnt, $orgs, $strs) = $ti->row($i);
    my @strs = split(" ", $strs);
    for (@strs) {
      my ($org, $id, $beg, $end, $srd) = split ":";
      print $fho join("\t", $cid, $org, $orgs, $id, $beg, $end, $srd)."\n";
    }
  }
  close $fho;
}
sub cluster_dedup {
  my ($fi, $fo) = @_;
  my $ti = readTable(-in => $fi, -header => 1);

  for my $i (0..$ti->lastRow) {
    my ($cid, $len, $cnt, $orgs, $strs) = $ti->row($i);
    my @strs = split(" ", $strs);
    my $ho;
    my @strs_n;
    for my $str (@strs) {
      my ($org, $id, $beg, $end, $srd) = split(":", $str);
      !exists $ho->{$org} || next;
      push @strs_n, $str;
      $ho->{$org} = 1;
    }
    $ti->setElm($i, "cnt", scalar(keys(%$ho)));
    $ti->setElm($i, "strs", join(" ", @strs_n));
  }
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $ti->tsv(1);
  close $fho;
}
sub cluster_filter {
  my ($fi, $fo, $minlen) = @_;
  my $ti = readTable(-in => $fi, -header => 1);

  my @idxs;
  for my $i (0..$ti->lastRow) {
    my ($cid, $len, $cnt, $orgs, $strs) = $ti->row($i);
    push @idxs, $i if $len < $minlen;
  }
  $ti->delRows(\@idxs);

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $ti->tsv(1);
  close $fho;
}
sub cluster_seq_prepare {
  my ($fi, $fo, $orgs_order) = @_;
  my $t = readTable(-in => $fi, -header => 1);

  my (@orgs, @ids, @begs, @ends);
  for my $i (0..$t->lastRow) {
    my ($cid, $len, $cnt, $orgs, $strs) = $t->row($i);
    my @strs = split(" ", $strs);
    my $hc;
    for my $str (@strs) {
      my ($org, $id, $beg, $end, $srd) = split(":", $str);
      !exists $hc->{$org} || die "$cid $strs has 2 $org\n";
      $hc->{$org} = [$id, $beg, $end, $srd];
    }

    my %ho = map {$_=>1} split(" ", $orgs);
    my @orgs_ary = grep {exists $ho{$_}} @$orgs_order;
    my $org = $orgs_ary[0];
    my ($id, $beg, $end, $srd) = @{$hc->{$org}};
    push @orgs, $org;
    push @ids,  $id;
    push @begs, $beg;
    push @ends, $end;
  }
  $t->addCol(\@orgs, 'org');
  $t->addCol(\@ids,  'id');
  $t->addCol(\@begs, 'beg');
  $t->addCol(\@ends, 'end');
  $t->sort("orgs", 1, 0, "id", 1, 0, "beg", 0, 0, "end", 0, 0); 

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $t->tsv(1);
  close $fho;
}
sub read_sv_tabix {
  my ($con, $id, $beg, $end) = @_;
  my $iter = $con->query($id, $beg, $end);
  my $flag = $iter->get();
  $flag || return ('chrX', 0);
  
  my @ary;
  while (my $line = $con->read($iter)) {
    my @ps = split("\t", $line);
    my ($tid, $tb, $te, $cid, $qid, $qb, $qe) = @ps;
    $beg > $tb || die "novseq[$id:$beg-$end] ovlp with aln $cid\[$tid:$tb-$te]\n";
    $end < $te || die "novseq[$id:$beg-$end] ovlp with aln $cid\[$tid:$tb-$te]\n";
    my $dist = min($beg - $tb, $te - $end);
    my $pos = int( ($qb + $qe) / 2 );
    push @ary, [$qid, $pos, $dist];
  }
  @ary = sort {$a->[2] <=> $b->[2]} @ary;
  my ($qid, $qpos) = @ary > 0 ? @{$ary[0]} : ('chrX', 0);
  return ($qid, $qpos);
}
sub seq_coord_anchor {
  my ($fi, $fo, $orgs) = @_;

  my $hg;
  for my $org (@$orgs) {
    my $fg = "/home/youngn/zhoup/Data/misc3/$org\_HM101/23_blat/41.9/sv.gz";
    my $con = Tabix->new(-data => $fg);
    $hg->{$org} = $con;
  }
  
  my $t = readTable(-in => $fi, -header => 1);
  my (@aids, @apos);
  for my $i (0..$t->lastRow) {
    my ($cid, $len, $cnt, $orgs2, $strs, $org, $id, $beg, $end) = $t->row($i);
    exists $hg->{$org} || die "no tabix for $org\n";
    my ($aid, $apos) = read_sv_tabix($hg->{$org}, $id, $beg, $end);
    push @aids, $aid;
    push @apos, $apos;
  }
  $t->addCol(\@aids,  'aid');
  $t->addCol(\@apos, 'apos');

  $t->sort("orgs", 1, 0, "aid", 1, 0, "apos", 0, 0, "id", 1, 0, "beg", 0, 0);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $t->tsv(1);
  close $fho;
}
sub cluster_seq {
  my ($fi, $fo, $fs) = @_;
  my $t = readTable(-in => $fi, -header => 1);

  my ($hs, $hd);
  my (@vids, @vbegs, @vends);
  for my $i (0..$t->lastRow) {
    my ($cid, $len, $cnt, $orgs, $strs, $org, $id, $beg, $end) = $t->row($i);
    if(!exists $hd->{$org}) {
      my $fd = "/home/youngn/zhoup/Data/genome/$org/11_genome.fas";
      my $db = Bio::DB::Fasta->new($fd);
      $hd->{$org} = $db;
    }
    my $vid = $orgs;
    $vid =~ s/ /\_/g;
    $hs->{$vid} ||= [0, ''];
    my ($pos, $pseq) = @{$hs->{$vid}};
    my ($vbeg, $vend) = ($pos + 1, $pos + $len);
    my $seq = $hd->{$org}->seq($id, $beg, $end);
    $pseq .= $seq . 'N' x 100;
    $pos += $len + 100;
    $hs->{$vid} = [$pos, $pseq];

    push @vids, $vid;
    push @vbegs, $vbeg;
    push @vends, $vend;
  }
  $t->addCol(\@vids,  'vid');
  $t->addCol(\@vbegs, 'vbeg');
  $t->addCol(\@vends, 'vend');

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $t->tsv(1);
  close $fho;
  
  my $seqH = Bio::SeqIO->new(-file=>">$fs", -format=>"fasta");
  for my $vid (sort(keys(%$hs))) {
    $seqH->write_seq(Bio::Seq->new(-id=>$vid, -seq=>$hs->{$vid}->[1]));
  }
  $seqH->close();
}

sub cluster_stat {
  print "total seq len:\n";
#  runCmd("fas2bed.pl -i 01.fas | bedlen.pl");
  
  my ($t, $len1, $len2);

  $t = readTable(-in => "13.clu", -header => 0);
  $len1 = sum( map {$t->elm($_, "col4") - $t->elm($_, "col3") + 1} (0..$t->lastRow) );
  printf "total clustered novseq len: %d\n", $len1;
  
  $t = readTable(-in => "21.cluster.tbl", -header => 1);
  $len1 = sum( map {$t->elm($_, "cnt") * $t->elm($_, "len")} (0..$t->lastRow) );
  $len2 = sum($t->col('len'));
  printf "21_cluster: all %d cluster %d\n", $len1, $len2;
  $t = readTable(-in => "21.coord.tbl", -header => 1);
  $len1 = sum( map {$t->elm($_, "end") - $t->elm($_, "beg") + 1} (0..$t->lastRow) );
  printf "21.coord: %d\n", $len1;

  $t = readTable(-in => "25.cluster.tbl", -header => 1);
  $len1 = sum( map {$t->elm($_, "cnt") * $t->elm($_, "len")} (0..$t->lastRow) );
  $len2 = sum($t->col('len'));
  printf "25_cluster: all %d cluster %d\n", $len1, $len2;
  $t = readTable(-in => "25.coord.tbl", -header => 1);
  $len1 = sum( map {$t->elm($_, "end") - $t->elm($_, "beg") + 1} (0..$t->lastRow) );
  printf "25.coord: %d\n", $len1;
 
  $t = readTable(-in => "27.cluster.tbl", -header => 1);
  $len1 = sum( map {$t->elm($_, "cnt") * $t->elm($_, "len")} (0..$t->lastRow) );
  $len2 = sum($t->col('len'));
  printf "27_cluster: all %d cluster %d\n", $len1, $len2;
  $t = readTable(-in => "27.coord.tbl", -header => 1);
  $len1 = sum( map {$t->elm($_, "end") - $t->elm($_, "beg") + 1} (0..$t->lastRow) );
  printf "27.coord: %d\n", $len1;
  
  $t = readTable(-in => "27.coord.tbl", -header => 1);
  my @orgstrs = ("HM340", "HM034", "HM056", "HM340 HM034", "HM340 HM056",
    "HM034 HM056", "HM340 HM034 HM056");
  for my $orgstr (@orgstrs) {
    my @orgs = split(" ", $orgstr);
    my $ho;
    for my $org (@orgs) {
      my $ts = $t->match_pattern_hash("\$_{'org'} eq '$org' && \$_{'orgs'} eq '$orgstr'");
      my $len = sum( map {$ts->elm($_, 'end') - $ts->elm($_, 'beg') + 1} (0..$ts->nofRow-1) );
      $ts = $ts->subTable([0..$ts->nofRow-1], [qw/id beg end/]);
      my $begs = $ts->delCol("beg");
      $begs = [ map {$_-1} @$begs ];
      $ts->addCol($begs, "beg", 1);
      open(my $fh, ">tmp1") or die "cannot write tmp1\n";
      print $fh $ts->tsv(0);
      close $fh;
      system("intersectBed -a \$genome/$org/51.bed/cds.bed -b tmp1 | bedSort stdin tmp2");
      my $lines = runCmd("mergeBed -i tmp2 | bedlen.pl", 2);
      my $lenc = $lines->[0];
      $ho->{$org} = [$len, $lenc];
      system("rm tmp1 tmp2");
    }
    my $str = join(" ", map {$_."[".$ho->{$_}->[1]."/".$ho->{$_}->[0]."]"} keys(%$ho));
    printf "$orgstr: $str\n";
  }
  
  $t = readTable(-in => "28.cluster.tbl", -header => 1);
  $len1 = sum( map {$t->elm($_, "cnt") * $t->elm($_, "len")} (0..$t->lastRow) );
  $len2 = sum($t->col('len'));
  printf "28_cluster: all %d cluster %d\n", $len1, $len2;
  
  $t = readTable(-in => "31.tbl", -header => 1);
  $len1 = sum( map {$t->elm($_, "cnt") * $t->elm($_, "len")} (0..$t->lastRow) );
  $len2 = sum($t->col('len'));
  my $len3 = 100 * $t->nofRow + $len2; 
  printf "31: all %d cluster %d seq[%d]\n", $len1, $len2, $len3;

  my $lines = runCmd("fas2bed.pl -i 33.pan4.fas | bedlen.pl", 2);
  printf "33.pan4: %d\n", $lines->[0];
}

sub clu_seq_blast {
  my ($fi, $fs, $fo, $minlen) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  my $db = Bio::DB::Fasta->new($fs);
  my $seqHO = Bio::SeqIO->new(-file => ">$fo", -format=>'fasta');
  for my $i (0..$t->nofRow-1) {
    my ($cid, $len, $org, $str_org, $str_cnt, $str_loc) = $t->row($i);
    next if $len < $minlen;
    my @strs = split(" ", $str_loc);
    my $idx = first_index {/^\Q$org\E/} @strs;
    my $str = $strs[$idx];
    my ($id, $beg, $end, $srd) = split(":", $str);
    my $seq = $db->seq($id, $beg, $end);
    $seq = Bio::Seq->new(-seq=>$seq)->revcom->seq if $srd eq "-";
    my $nid = "$cid|$str";
    $seqHO->write_seq(Bio::Seq->new(-id=>$nid, -seq=>$seq));
  }
  $seqHO->close();
}
sub merge_blastnr {
  my ($fi, $fb, $fo) = @_;
  my $t = readTable(-in => $fi, -header => 1);

  my $tb = readTable(-in => $fb, -header => 1);
  my $hb;
  for my $i (0..$tb->nofRow-1) {
    my ($id, $cat) = ($tb->elm($i, "id"), $tb->elm($i, "cat"));
    my ($cid, $str) = split(/\|/, $id);
    !exists $hb->{$cid} || die "$cid appeared >1 times in $fb\n";
    $hb->{$cid} = $cat;
  }

  my @srcs;
  for my $i (0..$t->nofRow-1) {
    my $cid = $t->elm($i, "cid");
    if(exists($hb->{$cid})) {
      push @srcs, $hb->{$cid};
    } else {
      push @srcs, "unc";
    }
  }

  $t->addCol(\@srcs, "src", 2);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $t->tsv(1);
  close $fho;
}



