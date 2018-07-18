#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use File::Path qw/make_path remove_tree/;
use Data::Dumper;
use List::Util qw/min max sum/;
use Common;
use Bam;

my $dir = "$ENV{'misc2'}/rnaseq/acc4";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir $dir\n";

#gzip_fastq();
run_tophat();

sub gzip_fastq {
  -d "12_fastq_gz" || make_path("12_fastq_gz");
  my $t = readTable(-in=>"01.tbl", -header=>1);
  for my $i (0..$t->lastRow) {
    next if $i < 36;
    my ($sam, $dir, $id, $label, $note1, $note2) = $t->row($i);
    my $fa = "$dir/$id\_$label\_R1_001.fastq";
    my $fb = "$dir/$id\_$label\_R2_001.fastq";
    -s $fa || die "$fa is not there\n";
    -s $fb || die "$fb is not there\n";
    print "compressing $id\_$label\n";
    runCmd("gzip -c $fa > 12_fastq_gz/$id.1.fq.gz", 1); 
    runCmd("gzip -c $fb > 12_fastq_gz/$id.2.fq.gz", 1); 
  }
}
sub read_fastq_by_org {
  my ($fi) = @_;
  my $t = readTable(-in=>$fi, -header=>1);
  my $h;
  for my $i (0..$t->lastRow) {
    my ($org, $dir, $id) = $t->row($i);
    $h->{$org} ||= [];
    push @{$h->{$org}}, $id;
  }
  return $h;
}
sub read_rnaseq_mapping {
  my ($fi) = @_;
  my $t = readTable(-in=>$fi, -header=>1);
  my $h;
  for my $i (0..$t->lastRow) {
    my ($org, $orgr) = $t->row($i);
    !exists $h->{$org} || die "$org has >=2 mappings\n";
    $h->{$org} = $orgr;
  }
  return $h;
}
sub run_tophat {
  -d "22_tophat" || make_path("22_tophat");
  my $hi = read_fastq_by_org("01.tsv");
  my $hr = read_rnaseq_mapping("21.tbl");
  for my $org (sort(keys(%$hr))) {
    my $orgr = $hr->{$org};
    my $dir = "22_tophat/$orgr\_$org";
    if(check_bam("$dir/unmapped.bam")) {
      printf "%s RNA-Seq -> %s: done\n", $orgr, $org;
      next;
    }
    printf "%s RNA-Seq -> %s: working...\n", $orgr, $org;
    my @ids = @{$hi->{$orgr}};
    my @fns;
    for my $id (@ids) {
      my $fa = "12_fastq_gz/$id.1.fq.gz";
      my $fb = "12_fastq_gz/$id.2.fq.gz";
      -s $fa || die "$fa is not there\n";
      -s $fb || die "$fb is not there\n";
      push @fns, ($fa, $fb);
    }
    my $str_fn = join(" ", @fns);
    runCmd("tophat2 --num-threads 24 --mate-inner-dist 0 \\
      --min-intron-length 45 --max-intron-length 5000 \\
      --min-segment-intron 45 --max-segment-intron 5000 \\
      -o $dir \$data/db/bowtie2/$org \\
      $str_fn");
    runCmd("samtools index $dir/accepted_hits.bam");
  }
}
sub merge_bam_genome {
  my ($fi, $di, $do) = @_;
  my $t = readTable(-in=>$fi, -header=>1);
  -d $do || make_path($do);

  my $h;
  for my $i (0..$t->lastRow) {
    my ($sam, $id, $genome) = $t->row($i);
    $h->{$genome} ||= [];
    push @{$h->{$genome}}, $id;
  }

  for my $gm (keys(%$h)) {
    my @ids = @{$h->{$gm}};
    for my $id (@ids) {
      check_bam("$di/$id\_$gm/accepted_hits.bam") || 
        die "no $di/$id\_$gm/accepted_hits.bam\n";
    }
    my $str_rg = "\@RG\\tID:$gm\\tSM:$gm\\tLB:$gm\\tPL:ILLUMINA\\tPU:lane";
    my $str_in = join(" ", map {"$di/$_\_$gm/accepted_hits.bam"} @ids);
    runCmd("samtools cat $str_in -o $do/$gm.bam");
  }
}
  
    

