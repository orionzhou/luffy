#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.expand.pl - 

=head1 SYNOPSIS
  
  comp.expand.pl [-help]

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
use Tabix;
use Bio::DB::Fasta;
use Bio::DB::Sam;
use Data::Dumper;
use Common;
use Location;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $tgt = "HM101";
my $qry = "HM340.AC";

my $dir = "$ENV{'misc3'}/comp.expand/$qry";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $f_ref = "$ENV{'genome'}/$qry/11_genome.fas";
my $f_bam = "$ENV{'misc3'}/pacbio/HM340_HM340.AC/15.bam";
find_pacbio_support("11.exp.tbl", "12.pacbio.tbl", $f_bam, $f_ref);

my $f_ort = "../../comp.ortho/33.ortho.cat.tbl";
my $f_tseq = "$ENV{'genome'}/$tgt/51.fas";
my $f_qseq = "$ENV{'genome'}/$qry/51.fas";
#make_subgroup_aln($f_ort, "05.aln", $tgt, $qry, $f_tseq, $f_qseq);

sub find_pacbio_support {
  my ($fi, $fo, $fb, $fr) = @_;
  my $ti = readTable(-in => $fi, -header => 1);

  my $bam = Bio::DB::Sam->new(-bam => $fb, -fasta => $fr);

  my ($nn, $nd, $np) = (0) x 3;
  my (@bposs, @eposs, @n_read, @reads);
  for my $i (0..$ti->lastRow) {
    my ($chr, $beg, $end, $srd, $id, $aid, $abeg, $aend, $bpos, $epos) = 
      $ti->row($i);
    $nn ++;
    if($aid eq "") {
      push @n_read, 0;
      push @reads, '';
      next;
    }
    $nd ++;
   
    my @alns = $bam->get_features_by_location(-seq_id => $chr, 
      -start => $bpos, -end => $epos);
    my @rns;
    for my $a (@alns) {
      my ($rbeg, $rend) = ($a->start, $a->end);
      if($rbeg <= $bpos && $rend >= $epos) {
        push @rns, $a->query->name;
      }
    }
    if(@rns > 0) {
      push @n_read, scalar(@rns);
      push @reads, join(" ", @rns);
      $np ++;
    } else {
      push @n_read, 0;
      push @reads, '';
    }
  }

  $ti->addCol(\@n_read, 'n_read');
  $ti->addCol(\@reads, 'reads');
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $ti->tsv(1);
  close $fho;

  printf "%d inserted/gained genes\n", $nn;
  printf "%d tandem duplications (with ancestors)\n", $nd;
  printf "%d with PacBio support\n", $np;
}
sub make_subgroup_aln {
  my ($fi, $do, $tgt, $qry, $tfas, $qfas) = @_;
  my $ti = readTable(-in => $fi, -header => 1);
  -d $do || make_path($do);
  runCmd("rm -rf $do/*");

  my $tdb = Bio::DB::Fasta->new($tfas);
  my $qdb = Bio::DB::Fasta->new($qfas);

  my $t = $ti->match_pattern_hash('$_{"cat2"} eq "CRP"');
  my ($ht, $hq);
  for my $i (0..$t->lastRow) {
    my $fam = $t->elm($i, 'cat3');
    my ($tid, $qid) = map {$t->elm($i, $_)} ($tgt, $qry);
    $ht->{$fam} ||= [];
    $hq->{$fam} ||= [];
    push @{$ht->{$fam}}, $tid if $tid ne "" && $tid ne "NA";
    push @{$hq->{$fam}}, $qid if $qid ne "" && $qid ne "NA";
  }
  my @fams = sort grep {@{$ht->{$_}} > 0 && @{$hq->{$_}} > 0} 
    uniq($t->col("cat3"));
  for my $fam (@fams) {
    print $fam."\n";
    my @tids = @{$ht->{$fam}};
    my @qids = @{$hq->{$fam}};
    open(my $fho, ">$do/$fam.fas") or die "cannot write $fam.fas\n";
    for my $tid (@tids) {
      my $seq = $tdb->seq($tid);
      print $fho ">$tgt|$tid\n$seq\n";
    }
    for my $qid (@qids) {
      my $seq = $qdb->seq($qid);
      print $fho ">$qry|$qid\n$seq\n";
    }
    close $fho;
    runCmd("clustalo -i $do/$fam.fas -o $do/$fam.aln \\
      --outfmt=clu --force --full --full-iter");
    runCmd("clustalw2 -infile=$do/$fam.aln -tree \\
      -outputtree=phylip -clustering=NJ -kimura");
  }
}
__END__

