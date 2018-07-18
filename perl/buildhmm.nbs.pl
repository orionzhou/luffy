#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Bio::SeqIO;
use Bio::DB::Fasta;
use Common;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;

my $dir = "$ENV{'misc2'}/nbs/mt_40";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

##### align mt4.0 nbs-lrrs
# gtb2fas.pl -i 01.gtb -d $genome/HM101/11_genome.fas -o 01.fas
# clustalo -i 01.fas -o 11.aln --outfmt=clu \
#   --iter=3 --guidetree-out=11.guide.nwk \
#   --output-order=tree-order --force --full --threads=16
# clustalw2 -infile=11.aln -bootstrap=1000 -outorder=aligned \ 
#   -outputtree=phylip -bootlabels=node -clustering=nj -kimura
# trimal -in 11.aln -out 12.aln -gappyout

#runCmd("aln2phy.pl -i 05.aln -o 05.phy -l 40");
#runCmd("phyml -i 05.phy -d aa");
#runCmd("mv 05.phy_phyml_tree.txt 06.nwk");
#runCmd("rm 05.phy_phyml*");

##### map mt1.0 nbs-lrrs to mt4.0
#extract_mt10_nbs_seq();
# blat -prot ../mt_40/01.fas 01.fas 11.psl
# psl2gal.pl -i 11.psl -o - | \
#   galbest.pl -i - -p qry -o - | \
#   galfilter.pl -i - -o 13.flt.gal -c 0.9 -m 30 -p 0.9 

##### phy.gene.R buld tree, manually annotate clade

#### build aln for subgroups
build_subgroup_aln('01.fas', '20.clade.tbl');
runCmd("cat 24_hmm/*.hmm > 30.all.hmm");

sub build_subgroup_aln {
  my ($fs, $fi) = @_;
  -d "21_fas" || make_path("21_fas");
  -d "22_aln" || make_path("22_aln");
  -d "23_trim" || make_path("23_trim");
  -d "24_hmm" || make_path("24_hmm");
  runCmd("rm 24_hmm/*");
  my $db = Bio::DB::Fasta->new($fs);
  
  my $h;
  my $t = readTable(-in => $fi, -header => 1);
  for my $i (0..$t->lastRow) {
    my ($id, $oclade, $clade) = $t->row($i);
    $clade ne '' || next;
    my $seq = $db->seq($id);
    $h->{$clade} ||= [];
    push @{$h->{$clade}}, [$id, $seq];
  }
  for my $cla (sort(keys(%$h))) {
    $cla =~ /^[tc]nl/i || next;
    open(my $fhc, ">21_fas/$cla.fas") || die "cannot write 21_fas/$cla.fas\n";
    my @ary = @{$h->{$cla}};
    for (@ary) {
      my ($id, $seq) = @$_;
      $id =~ s/Medtr//i;
      my $nid = "$cla\_$id";
      print $fhc ">$nid\n$seq\n";
    }
    close $fhc;
    runCmd("clustalo -i 21_fas/$cla.fas -o 22_aln/$cla.aln \\
      --outfmt=clu --output-order=tree-order --force --full --threads=16");
    runCmd("trimal -automated1 -in 22_aln/$cla.aln -out 23_trim/$cla.aln");
    runCmd("hmmbuild 24_hmm/$cla.hmm 23_trim/$cla.aln"); 
  }
}

sub extract_mt10_nbs_seq {
  my $dir = "$ENV{'misc2'}/nbs/mt_10";
  my $fi = "$dir/00.raw.fas";
  my $fo = "$dir/01.fas";
  my $shi = Bio::SeqIO->new(-file => "<$fi", -format => 'fasta');
  my $sho = Bio::SeqIO->new(-file => ">$fo", -format => 'fasta');
  while( my $seqo = $shi->next_seq() ) {
    my $id = $seqo->id;
    if($id =~ /^(Mt[0-9a-z]+)/i) {
      my $nid = $1;
      $sho->write_seq(Bio::Seq->new(-id => $nid, -seq => $seqo->seq));
    } else {
      next;
    }
  }
  $shi->close();
  $sho->close();
}

sub extract_seq_by_pfam {
  my ($fo, $orgs, $fam) = @_;
  open(my $fho, ">$fo") or die "cannot write $fo\n";
 
  my $di = $ENV{'genome'};
  for my $org (@$orgs) {
    my $fd = "$di/$org/21.fas"; 
    my $fp = "$di/$org/23.pfam.tsv"; 
    my $tp = readTable(-in => $fp, -header => 0);
    my $db = Bio::DB::Fasta->new($fd);
  
    for my $i (0..$tp->nofRow-1) {
      my ($id, $md5, $len, $prog, $pid, $desc, $beg, $end, $e) = $tp->row($i);
      if($pid eq $fam) {
        my $nid = "$org-$id-$beg-$end";
        my $seq = $db->seq($id, $beg, $end);
        print $fho ">$nid\n$seq\n";
      }
    }
  }
  close $fho;
}



