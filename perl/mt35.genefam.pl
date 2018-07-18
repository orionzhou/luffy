#!/usr/bin/perl -w
use strict;
use Common;
use Gtb;
use Data::Dumper;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $dir = "$ENV{'misc2'}/genefam";
my $f01 = "$dir/01.gtb";
my $f11 = "$dir/11_crp.gtb";
my $f12 = "$dir/12_nbs.gtb";
my $f21 = "$dir/21_clean.gtb";
setDiffGtb(-in=>$f01, -crp=>$f11, -nbs=>$f12, -out=>$f21);
my $f31 = "$dir/31_merged.gtb";
#runCmd("awk '{ if( (NR==1) || (NR>1 && \$1!=\"id\") ) print }' $f21 $f11 $f12 > $f31";

my $f41 = "$dir/41_genefam.tbl";
my $f35 = "$dir/35";
#gtb_split_cat($f31, $f41, $f35);
sub gtb_split_cat {
  my ($fi, $f_cat, $fo) = @_;
  my $ti = readTable(-in=>$fi, -header=>1);
  my $tc = readTable(-in=>$f_cat, -header=>1);
  my @cats = uniq($tc->col("cat"));
  my $h = { map {$tc->elm($_, "id") => $tc->elm($_, "cat")} (0..$tc->nofRow-1) };
  for my $cat (@cats) {
    my $f_gtb = "$fo\_$cat.gtb";
    open(FH, ">$f_gtb");
    print FH join("\t", $ti->header)."\n";
    for my $i (0..$ti->nofRow-1) {
      my $id = $ti->elm($i, "id");
      die "$id not in $f_cat\n" unless exists $h->{$id};
      next unless $h->{$id} eq $cat;
      print FH join("\t", $ti->row($i))."\n";
    }
    close FH;
    my $f_gff = "$fo\_$cat.gff";
    gtb2Gff($f_gtb, $f_gff, 2);
  }
}
sub setDiffGtb {
  my ($fi, $f_crp, $f_nbs, $fo) = rearrange(['in', 'crp', 'nbs', 'out'], @_);

  my @ids_crp;
  my $t_crp = readTable(-in=>$f_crp, -header=>1);
  
  my $f_crp_info = "$ENV{'misc2'}/crp/05_models/01_mt35_v5/01_ovlp_gene.tbl";
  my $tt = readTable(-in=>$f_crp_info, -header=>1);
  my $h = { map {$tt->elm($_, "id") => $tt->elm($_, "gene")} (0..$tt->nofRow-1) };

  for my $i (0..$t_crp->nofRow-1) {
    my ($id_crp, $tag) = map {$t_crp->elm($i, $_)} qw/parent conf/;
    $id_crp =~ /CRP\d+\_(\d+)/i;
    my $id_crp_num = int($1);
    next unless exists $h->{$id_crp_num};
    push @ids_crp, $h->{$id_crp_num};
  }
  my $h_crp = { map {$_=>1} @ids_crp };
  printf "%d IDs from CRP\n", scalar(keys(%$h_crp));

  my $t_nbs = readTable(-in=>$f_nbs, -header=>1);
  my $h_nbs = { map {$_=>1} $t_nbs->col("id") };
  printf "%d IDs from NBS-LRR\n", scalar(keys(%$h_nbs));

  my $t = readTable(-in=>$fi, -header=>1);
  my @idxs;
  for my $i (0..$t->nofRow-1) {
    my $id = $t->elm($i, "id");
    if(exists $h_crp->{$id} || exists $h_nbs->{$id}) {
      push @idxs, $i;
    }
  }
  $t->delRows(\@idxs);
  printf "%d IDs deleted\n", scalar(@idxs);
  open(FH, ">$fo");
  print FH $t->csv(1, {delimiter=>"\t"});
  close FH;
}


