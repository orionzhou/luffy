#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  mt.augus.pl - run augustus in parallel on Medicago genome assembly

=head1 SYNOPSIS
  
  mt.augus.pl [-help] [-org organism]

  Options:
    -h (--help)   brief help message
    -g (--org)    genmome ID (organism) to process

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Location;
use Data::Dumper;
use Gtb;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($org, $ncpu) = ('', 24);
my $help_flag;
#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "org|g=s"  => \$org,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$org;

my $dir = "$ENV{'genome'}/$org/augustus";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir $dir\n";

my $fg = "../11_genome.fas";
-s $fg || die "$fg is not there\n";

get_hints_ortholog();
get_hints_rnaseq();
run_aug();
postprocess_aug();
pipe_pfam();

sub get_hints_ortholog {
  my $f_gtb = "$ENV{'genome'}/HM101/51.gtb";
  my $f_ref = "$ENV{'genome'}/HM101/11_genome.fas";
  my $f_gax = "$ENV{'misc3'}/$org\_HM101/23_blat/31.9/gax.gz";
  my $f_snp = "$ENV{'misc3'}/$org\_HM101/23_blat/31.9/snp.gz";
  runCmd("liftover.gtb2hint.pl -i $f_gtb -r $f_ref -x $f_gax -s $f_snp -o 05.hm101.gff");
}
sub get_hints_rnaseq {
  my $dirin = "$ENV{'misc2'}/rnaseq/acc4";
  
  my $t = readTable(-in => "$dirin/21.tbl", -header => 1);
  my $h;
  for my $i (0..$t->lastRow) {
    my ($org, $orgr) = $t->row($i);
    !exists $h->{$org} || die "$org has >=2 mappings\n";
    $h->{$org} = $orgr;
  }
  exists $h->{$org} || die "no RNA-seq for $org\n";
  my $orgr = $h->{$org};
  
  my $bam_in = "$dirin/22_tophat/$orgr\_$org/accepted_hits.bam";
  -s $bam_in || die "$bam_in is not there\n";
  runCmd("bamtools filter -isPrimaryAlignment 1 -isMapped 1 -isMateMapped 1 -in $bam_in -out 11.f.bam");
  runCmd("samtools view -H 11.f.bam > 12.header.txt");

  runCmd("samtools sort 11.f.bam -o 14.sf.bam");
  runCmd("bam2hints --intronsonly --in=14.sf.bam --out=15.rnaseq.gff");
#  runCmd("bam2wig $bam_in | \$soft/augustus/scripts/wig2hints.pl > 15.ep.gff");

#  runCmd("rm 11.f.bam 12.header.txt 14.sf.bam");
}
sub run_aug {
  -s "05.hm101.gff" || die "no 05.hm101.gff";
  -s "15.rnaseq.gff" || die "no 15.rnaseq.gff";
  runCmd("cat 05.hm101.gff 15.rnaseq.gff > 20.hints.gff");
  runCmd("sort -k1,1 -k4,4n 20.hints.gff -o 20.hints.gff");

  my $dig = getDigits($ncpu);
  runCmd("ln -sf $fg 21.fas");
  runCmd("pyfasta split -n $ncpu 21.fas");
  my $cfg_aug = "$ENV{'AUGUSTUS_CONFIG_PATH'}/extrinsic/extrinsic.mt.cfg";
  my $end = $ncpu - 1;
  runCmd("seq 0 $end | xargs -i printf \"%0".$dig."d\\n\" {} | \\
    parallel -j $ncpu augustus --species=medicago \\
    --extrinsicCfgFile=$cfg_aug \\
    --alternatives-from-evidence=true \\
    --allow_hinted_splicesites=atac \\
    --introns=on --genemodel=partial \\
    --strand=both --gff3=on \\
    --hintsfile=20.hints.gff \\
    --outfile=21.{}.gff 21.{}.fas");
  runCmd("cat 21.*.gff | join_aug_pred.pl > 21.gff");
  runCmd("rm 21.fas.* 21.*.fas 21.*.gff");
}
sub postprocess_aug {
  runCmd("gff.augus.pl -i 21.gff -o - | gff2gtb.pl -i - -o 22.gtb");
  runCmd("gtb.rmutr.pl -i 22.gtb -o 23.gtb");
  runCmd("gtb.dedup.pl -i 23.gtb -o 25.dedup.gtb");
  runCmd("gtb2gtbx.pl -i 25.dedup.gtb -d $fg -o 31.gtbx");
  runCmd("cut -f1-18 31.gtbx > 31.gtb");
  runCmd("gtbx2fas.pl -i 31.gtbx -o 31.fas");
  runCmd("gtb2gff.pl -i 31.gtb -o 31.gff");
}
sub pipe_pfam {  
##runCmd("interproscan.sh -appl PfamA -dp -i 31.fas -f tsv -o 33.pfam.tsv");
##runCmd("pfam2tbl.pl -i 33.pfam.tsv -o 34.pfam.tbl -e 1 -l 10");
  #-s "33.txt" && runCmd("cp 33.txt 34.tbl.1.txt");
  runCmd("pfam.scan.pl -i 31.fas -o 34.tbl");
  runCmd("gtb.addpfam.pl -i 31.gtb -p 34.tbl -o 41.gtb"); 
  gtb2bed("41.gtb", "41.bed");
  runCmd("intersectBed -wao -a 41.bed -b ../12.rm.bed > 42.ovlp.bed");
  gtb_add_rm("41.gtb", "42.ovlp.bed", "43.gtb");
}
sub pipe_blastnr {
#  runCmd("awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {if(NR>1 && \$16 == \"Unknown\") print \$1}' 41.gtb > 42.unk.txt");
#  runCmd("seqret.pl -d 31.fas -b 42.unk.txt -o 42.unk.fas");
  runCmd("pro.blastnr.py 42.unk.fas 44");
}
sub gtb2bed {
  my ($fi, $fo) = @_;
  open(my $fhi, $fi) || die "cannot read file $fi\n";
  open(my $fho, ">$fo") || die "cannot write $fo\n";
  while(<$fhi>) {
    chomp;
    /(^id)|(^\#)|(^\s*$)/i && next;
    my $ps = [ split("\t", $_, -1) ];
    @$ps >= 18 || die "not 19 fileds:\n$_\n";
    my ($id, $par, $chr, $beg, $end, $srd, 
      $locES, $locIS, $locCS, $loc5S, $loc3S, $phase, 
      $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
    print $fho join("\t", $chr, $beg-1, $end, $id)."\n";
  }
  close $fhi;
  close $fho;
}
sub gtb_add_rm {
  my ($fi, $fb, $fo) = @_;
  open(my $fhb, $fb) || die "cannot read file $fb\n";
  my $h = {};
  while(<$fhb>) {
    chomp;
    my ($c1, $b1, $e1, $id, $c2, $b2, $e2) = split "\t";
    $h->{$id} ||= 0;
    $h->{$id} += $e2 - $b2;
  }
  close $fhb;
  
  open(my $fhi, $fi) || die "cannot read file $fi\n";
  open(my $fho, ">$fo") || die "cannot write $fo\n";
  print $fho join("\t", @HEAD_GTB)."\n";
  my $cnt = 0;
  my $h2 = {};
  while(<$fhi>) {
    chomp;
    /(^id)|(^\#)|(^\s*$)/i && next;
    my $ps = [ split("\t", $_, -1) ];
    @$ps >= 18 || die "not 19 fileds:\n$_\n";
    my ($id, $par, $chr, $beg, $end, $srd, 
      $locES, $locIS, $locCS, $loc5S, $loc3S, $phase, 
      $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
    my $leng = $end - $beg + 1;
    my $lenr = exists $h->{$id} ? $h->{$id} : 0;
    if($lenr / $leng >= 0.6 && $cat2 ne 'TE') {
      $ps->[15] = 'TE';
      $ps->[16] = 'RepeatMasker';
      $cnt ++;
      $h2->{$cat2} ||= 0;
      $h2->{$cat2} ++;
    }
    print $fho join("\t", @$ps)."\n";
  }
  close $fhi;
  close $fho;
  print "$cnt new TEs added\n";
  for my $fam (keys(%$h2)) {
    print "$fam\t$h2->{$fam}\n";
  }
}


__END__

