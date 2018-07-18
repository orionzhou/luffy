#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  pipe.vntcall.pl - call SNPs and coverage from Bam and Vcf files

=head1 SYNOPSIS
  
  pipe.vntcall.pl [-help] 

  Options:
    -h (--help)   brief help message

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common; 
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use Cwd qw/abs_path/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my @orgs = qw/
  HM058 HM125 HM056 HM129 HM060
  HM095 HM185 HM034 HM004 HM050 
  HM023 HM010 HM022 HM324 HM340
/;
#@orgs = qw/
#  HM058 HM125 HM056 HM129 HM060
#  HM095 HM185 HM034 HM004 HM050 
#  HM023 HM010/;

my $dir = "$ENV{'misc3'}/hapmap_mt40/12_ncgr";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $fs = "$ENV{'genome'}/HM101/15.sizes";
my $fv = "$dir/../30_vnt/acc319.vcf.gz";
-s $fv || die "cannot find $fv\n";

-d "35_cov" || make_path("35_cov");
-d "36_abcov" || make_path("36_abcov");
-d "44_snp" || make_path("44_snp");

#call_vcf();
#call_vnt_org();
#call_indel();
call_cvg();

sub call_vcf {
  my @norgs;
  for (@orgs) { 
    $_ =~ s/^HM0((17)|(20)|(22))$/HM0$1-I/;
    push @norgs, $_;
  }
  my $orgstr = join(",", @norgs);
  my $chrstr = join(",", map {"chr".$_} (1..8));
  runCmd("bcftools view -r $chrstr \\
    --exclude-uncalled -g hom --min-ac 2:alt1 --trim-alt-alleles \\
    --output-type v --types snps,indels,mnps \\
    -s $orgstr $fv -o 41.vcf");
  runCmd("bcftools norm -m +snps --check-ref w -f \$genome/HM101/11_genome.fas -O v 41.vcf -o 42.norm.vcf");
  runCmd("bcftools view -v snps -O v 42.norm.vcf | vcf.addref.pl -i - -o 43.snp.vcf");
  runCmd("vcf.stat.py 43.snp.vcf 45.stat.tbl");
}
sub call_vnt_org {
  for my $org (@orgs) {
    my $orgi = $org;
    $orgi =~ s/^HM0((17)|(20)|(22))$/HM0$1-I/;
    runCmd("bcftools view \\
      --exclude-uncalled --min-ac 1:alt1 --trim-alt-alleles \\
      -s $orgi -O v 43.snp.vcf -o 44_snp/$org.vcf");
    runCmd("vcf2tbl.py 44_snp/$org.vcf 44_snp/$org.tbl");
    runCmd("bgzip -c 44_snp/$org.tbl > 44_snp/$org.tbl.gz");
    runCmd("tabix -s 1 -b 2 -e 2 -S 1 44_snp/$org.tbl.gz");
  }
}
sub call_indel {
  my @norgs;
  for (@orgs) { 
    $_ =~ s/^HM0((17)|(20)|(22))$/HM0$1-I/;
    push @norgs, $_;
  }
  my $orgstr = join(",", @norgs);
  my $chrstr = join(",", map {"chr".$_} (1..8));
  runCmd("bcftools view -v indels -O v 42.norm.vcf | vcf.addref.pl -i - -o 51.indel.vcf");
  runCmd("vcf.stat.py 51.indel.vcf 52.stat.tbl");
}
sub call_cvg {
  for my $org (@orgs) {
    my $orgi = $org;
    $orgi =~ s/^HM0((17)|(20)|(22))$/HM0$1-I/;
    
    my $fb = "01_raw/$orgi.recalibrated.bam";
    -s $fb || die "cannot find $fb\n";
   
  runCmd("bamtools filter -in $fb -mapQuality \">=20\" \\
    -isDuplicate \"false\" | samtools depth /dev/stdin | \\
    cov2bed.pl --chr-only -o 35_cov/$org.raw.bed");
  runCmd("bedGraphToBigWig 35_cov/$org.raw.bed $fs 35_cov/$org.bw");
  runCmd("mergeBed -i 35_cov/$org.raw.bed > 35_cov/$org.bed");

#  runCmd("bamFilterAb -i $fb -o 36_abcov/$org.bam");
#  runCmd("bamtools filter -in 36_abcov/$org.bam -mapQuality \">=20\" \\
#    -isDuplicate \"false\" | samtools depth /dev/stdin | \\
#    cov2bed.pl -o 36_abcov/$org.bed");
#  runCmd("bedGraphToBigWig 36_abcov/$org.bed $fs 36_abcov/$org.bw");
#  runCmd("rm 36_abcov/$org.bed");
  }
}

__END__
sub collect_rg {
  my ($ids, $d01, $f02, $d09) = @_;
  my $t = Data::Table->new([], [qw/sm rg lb id_lb id_run pl rl pi fq1 fq2 ds fn/]);
  for my $id (@$ids) {
    my $fi = "$d01/$id.bam";
    die "$fi is not there\n" unless -s $fi;
    open(PP, "samtools view -H $fi | grep '\@RG' |") or die "error piping\n";
    my $hl = {};
    while( <PP> ) {
      chomp;
      my @ps = split "\t";
      next unless $ps[0] eq '@RG';
      my ($fn, $pi, $ds) = ("") x 3;
      for my $ele (@ps[1..$#ps]) {
        my ($tag, $value) = split(":", $ele); 
        $fn = $value if $tag eq "ID";
        $pi = $value if $tag eq "PI";
        $ds = $value if $tag eq "DS";
      }
      my ($id_lb, $id_run) = (1, 1);
      if(exists($hl->{$pi})) {
        ($id_lb, $id_run) = @{$hl->{$pi}};
        $hl->{$pi}->[1] = ++$id_run;
      } else {
        $id_lb = scalar(keys(%$hl)) + 1;
        $hl->{$pi} = [$id_lb, $id_run];
      }
      my $lb = sprintf "%s_%02d", $id, $id_lb;
      my $rg = sprintf "%s_%02d_%02d", $id, $id_lb, $id_run;
      my ($fq1, $fq2) = map {"$d09/$rg.$_.fq.gz"} (1,2);
      $t->addRow([$id, $rg, $lb, $id_lb, $id_run, 'ILLUMINA', 90, $pi, $fq1, $fq2, $ds, $fn]);
    }
    printf "  collecting \@RG tags for %s...\r", $id;
  }
  print "\n";
  $t->sort('rg', 1, 0);
  open(FH, ">$f02");
  print FH $t->tsv(1);
  close FH;
}
sub pipe_ncgr {
  my ($dir, $beg, $end) = @_;
  my $d01 = "$dir/01_raw";
  my $f02 = "$dir/02_readgroups.tbl";
  my $d03 = "$dir/03_readname_sorted";
  my $d09 = "$dir/09_fastq";
  for my $i ($beg..$end) {
    my $id = sprintf "HM%03d", $i;
    my $f01 = "$d01/$id.bam";
#  runCmd("samtools sort -m 5000000000 -n $f01 $d03/$id");
    my $f03 = "$d03/$id.bam";
    runCmd("bam2Fastq -i $f03 -o $d09 -s $id -m $f02", 1);
#    runCmd("samtools sort -m 5000000000 -n $d01/$id.bam $d03/$id", 1);
#    runCmd("bamPreprocess -i $f03 -o $d04/$id", 1);
#    runCmd("samtools sort -m 5000000000 $d04/$id $d05/$id", 1);
#    runCmd("samtools index $d05/$id.bam", 1);
  }
}

sub merge_fq_list {
  my ($fi, $f_add, $fo) = @_;
  my $ti = readTable(-in=>$fi, -header=>1);
  my $ta = readTable(-in=>$f_add, -header=>1);
  my $h;
  for my $i (0..$ti->nofRow-1) {
    my ($sm, $rg, $lb, $id_lb, $id_run, $pl, $rl, $pi, $fq1, $fq2, $ds, $fn) = $ti->row($i);
    die "[$rg/fastq1] $fq1 is not there\n" unless -s $fq1;
    die "[$rg/fastq2] $fq2 is not there\n" unless -s $fq2;
    $h->{$sm}->{$pi} ||= [$id_lb, []];
    push @{$h->{$sm}->{$pi}->[1]}, $id_run;
  }
  for my $i (0..$ta->nofRow-1) {
    my ($sm, $rg, $lb, $id_lb, $id_run, $pl, $rl, $pi, $fq1, $fq2, $ds, $fn) = $ta->row($i);
    die "[$rg/fastq1] $fq1 is not there\n" unless -s $fq1;
    die "[$rg/fastq2] $fq2 is not there\n" unless -s $fq2;
    if(exists($h->{$sm}->{$pi})) {
      $id_lb = $h->{$sm}->{$pi}->[0];
      my $ids_run = $h->{$sm}->{$pi}->[1];
      $id_run = max(@$ids_run) + 1;
      push @$ids_run, $id_run;
    } else {
      my @ids_lb = map {$_->[0]} values(%{$h->{$sm}});
      $id_lb = max(@ids_lb) + 1;
      $id_run = 1;
      $h->{$sm}->{$pi} = [$id_lb, [$id_run]];
    }
    $pl ||= "ILLUMINA";
    $lb = sprintf "%s_%02d", $sm, $id_lb;
    $rg = sprintf "%s_%02d_%02d", $sm, $id_lb, $id_run;
    $ti->addRow([$sm, $rg, $lb, $id_lb, $id_run, $pl, $rl, $pi, $fq1, $fq2, $ds, $fn]);
  }
  $ti->sort("rg", 1, 0);
  $ti->addCol([1..$ti->nofRow], "idx", 0);
  open(FH, ">$fo");
  print FH $ti->tsv(1);
  close FH;
}
sub process_fastq {
  my ($d10, $beg, $end) = @_;
  for ($beg..$end) {
    my $acc = sprintf "HM%03d", $_;
    my $fi1 = "$d10/$acc/s_3_1_sequence.txt.gz";
    my $fi2 = "$d10/$acc/s_3_2_sequence.txt.gz";
    open(my $fhi1, "zcat $fi1 |") or die "cannot open $fi1\n";
    open(my $fhi2, "zcat $fi2 |") or die "cannot open $fi2\n";
    my $i = 0;
    my $h = {};
    while( !eof($fhi1) && !eof($fhi2) ) {
      my (@ls1, @ls2);
      for my $i (0..3) {
        push @ls1, substr(readline($fhi1), 0, -1);
        push @ls2, substr(readline($fhi2), 0, -1);
      }
      $ls1[0] =~ /^\@([\w\:\#]+)\/([12]) run=([\w\_]+)/;
      my ($id, $pair, $run) = ($1, $2, $3);
      $ls2[0] =~ /^\@([\w\:\#]+)\/([12]) run=([\w\_]+)/;
      die "error matching reads \n $id $1 $pair $2 $run $3\n" unless $id eq $1 && $run eq $3 && $pair == 1 && $2 == 2;
      $ls1[2] =~ /^\+([\w\:\#]+)\/([12])/;
      die "error matching reads\n" unless $id eq $1 && $2 == 1; 
      $ls2[2] =~ /^\+([\w\:\#]+)\/([12])/;
      die "error matching reads\n" unless $id eq $1 && $2 == 2; 

      if(exists($h->{$run})) {
        my ($fho1, $fho2) = @{$h->{$run}};
        print $fho1 join("\n", "\@$id/1", $ls1[1], "+", $ls1[3])."\n";
        print $fho2 join("\n", "\@$id/1", $ls2[1], "+", $ls2[3])."\n";
      } else {
        $i ++;
        my $fo1 = sprintf "$d10/$acc\_%02d.1.fq", 20+$i;
        my $fo2 = sprintf "$d10/$acc\_%02d.2.fq", 20+$i;
        my $fho1 = new IO::File $fo1, "w";
        my $fho2 = new IO::File $fo2, "w";
        $h->{$run} = [$fho1, $fho2];
      }
    }
    for my $j (1..$i) {
      my $fo1 = sprintf "$d10/$acc\_%02d.1.fq", 20+$j;
      my $fo2 = sprintf "$d10/$acc\_%02d.2.fq", 20+$j;
      runCmd("gzip $fo1");
      runCmd("gzip $fo2");
    }
  }
}
sub run_mpileup {
  my ($dir_bam, $d01, $d02, $ids) = @_;
  my @chrs = map {"chr$_"} (1..8);
#  @chrs = map {"chr$_"} (5);
  for my $id (@$ids) {
    print "working on $id\n";
    my $f00 = file($dir_bam, "$id.bam");
    my $f01 = file($d01, "$id.bcf");
    my $f02 = file($d02, "$id.bcf");

    my @f1_chrs;
    my @f2_chrs;
    for my $chr (@chrs) {
      my $f1_chr = file($d01, "$id\_$chr.bcf");
      my $f2_chr = file($d02, "$id\_$chr.bcf");
      runCmd2("samtools mpileup -ugf $f_genome -r $chr $f00 | bcftools view -bcg - > $f1_chr");
      runCmd2("samtools mpileup -q1 -ugf $f_genome -r $chr $f00 | bcftools view -bcg - > $f2_chr");
      push @f1_chrs, $f1_chr;
      push @f2_chrs, $f2_chr;
    }
    runCmd2("bcftools cat ".join(" ", @f1_chrs)." > $f01");
    runCmd2("bcftools index $f01");
    runCmd2("bcftools cat ".join(" ", @f2_chrs)." > $f02");
    runCmd2("bcftools index $f02");
    runCmd2("rm ".join(" ", @f1_chrs));
    runCmd2("rm ".join(" ", @f2_chrs));
  }
}



