#!/usr/bin/perl
use strict;
use Common; 
use Seq;
use Data::Dumper;

my $dir = "$DIR_Misc1/stats_gene";
my $f_seq = "$DIR_Genome/mt_35/41_genome.fa";
my $f_gtb = "$DIR_Genome/mt_35/10_model_Mt3.5v5/62_phase_fixed.gtb";
my $f01 = "$dir/01_gc.tbl";
get_gene_gc($f_gtb, $f_seq, $f01);

sub get_gene_gc {
  my ($f_gtb, $f_seq, $fo) = @_;
  my $tg = readTable(-in=>$f_gtb, -header=>1);
  open(FH, ">$fo");
  print FH join("\t", qw/id length a t g c n pct_gc/)."\n";
  for my $i (0..$tg->nofRow-1) {
    my ($id, $chr, $srd, $locCS) = map {$tg->elm($i, $_)} qw/id chr strand locC/;
    my $locC = locStr2Ary($locCS);
    my $seq = seqRet($locC, $chr, $srd, $f_seq);
    my $len = locAryLen($locC);
    die "$id not $len bp: ".length($seq)."\n" unless length($seq) == $len;
    my ($a, $t, $g, $c, $n) = (0) x 5;
    for my $i (0..$len-1) {
      my $ch = substr($seq, $i, 1);
      if($ch =~ /a/i) {
        $a ++;
      } elsif($ch =~ /t/i) {
        $t ++;
      } elsif($ch =~ /c/i) {
        $c ++;
      } elsif($ch =~ /g/i) {
        $g ++;
      } else {
        $n ++;
      }
    }
    print FH join("\t", $id, $len, $a, $t, $g, $c, $n, sprintf("%.03f", ($c+$g)/($len-$n)))."\n";
    printf "  %6d out of %6d done...\r", $i+1, $tg->nofRow;
  }
  print "\n";
  close FH;
}


__END__

