#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  panseq.pl - construct pan-genome fasta

=head1 SYNOPSIS
  
  panseq.pl [-help]

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
use Bio::AlignIO;
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

my $dir = "$ENV{'misc3'}/comp.panseq";
chdir $dir || die "cannot chdir to $dir\n";

my @orgs = qw/
  HM058 HM125 HM056 HM129 HM060
  HM095 HM185 HM034 HM004 HM050 
  HM023 HM010 HM022 HM324 HM340
/;

if($stat_flag) {exit;}

#write_file_list("01_files", "02.file.list", \@orgs);

my $dir_tmp = "/lustre/zhoup/tmp_paramugsy";
my $f_tmpl = "\$soft/paramugsy/pm_qsub_template.sh";
#runCmd("paramugsy local -cores 16 -seq-list 01.filelist.txt \\
#  -out-maf 21.maf -tmp-dir $dir_tmp -template-file $f_tmpl");

#maf2tbl('21.maf', '22.tbl');
restore_coordinate("31.refined.tbl", "32.global.tbl");
sub write_file_list {
  my ($do, $fl, $orgs) = @_;
  -d $do || make_path($do);
  open(my $fhl, ">$fl") or die "cannot write $fl\n";
  for my $org (@$orgs) {
    my $fn = "$ENV{'misc3'}/$org\_HM101/41_novseq/21.fas";
    runCmd("ln -sf $fn $do/$org.fas");
    print $fhl "$dir/$do/$org.fas\n";
  }
  close $fhl;
}
sub maf2tbl {
  my ($fi, $fo) = @_;
  my $ai = Bio::AlignIO->new(-file => $fi, -format => 'maf');
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/cid alen org chr beg end srd/)."\n";
  my $i = 1;
  while(my $aln = $ai->next_aln()) {
    my $n = $aln->num_sequences;
    my $alen = $aln->length();
    for my $seq ($aln->each_seq) {  
      my ($id, $beg, $end, $srd) = ($seq->id, $seq->start, $seq->end,
        $seq->strand);
      my ($org, $nid) = split(/\./, $id);
      $srd = $srd == -1 ? "-" : "+";
      if($srd eq "-" && $end > $beg) {
        my ($pid, $pbeg, $pend) = split("_", $nid);
        my $len = $pend - $pbeg + 1;
        ($beg, $end) = ($len - $end + 1, $len - $beg + 1);
      }
      print $fho join("\t", $i, $alen, $org, $nid, $beg, $end, $srd)."\n";
    }
    $i ++;
  }
  close $fho;
}
sub restore_coordinate {
  my ($fi, $fo) = @_;
  my $ti = readTable(-in => $fi, -header => 1);
  for my $i (0..$ti->lastRow) {
    my ($cid, $alen, $org, $sid, $rb, $re, $srd) = $ti->row($i);
    my ($chr, $bb, $be) = split("_", $sid);
#    my $beg = $srd eq "-" ? $be - $re + 1 : $bb + $rb - 1;
#    my $end = $srd eq "-" ? $be - $rb + 1 : $bb + $re - 1;
    my ($beg, $end) = ($bb + $rb - 1, $bb + $re - 1);
    $beg < $bb && die "error: $sid:$rb-$re => $chr:$beg:$end\n";
    $end > $be && die "error: $sid:$rb-$re => $chr:$beg:$end\n";
    $ti->setElm($i, "chr", $chr);
    $ti->setElm($i, "beg", $beg);
    $ti->setElm($i, "end", $end);
  }
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $ti->tsv(1);
  close $fho;
}

