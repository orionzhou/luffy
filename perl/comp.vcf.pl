#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.vcf.pl - merge VCF forat

=head1 SYNOPSIS
  
  comp.vcf.pl [-help] 

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
use Data::Dumper;
use Common;
use Location;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;
#pod2usage(2) if !$qry || !$tgt;

my $dir = "$ENV{'data'}/misc3/comp.vnt/acc4";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $tgt = "HM101";
my @qrys = qw/
  HM058 HM125 HM056 HM129 HM060
  HM095 HM185 HM034 HM004 HM050 
  HM023 HM010/;
@qrys = qw/HM056 HM034 HM340/;

#merge_vcfs(\@qrys);
#vcf2tbl("01.raw.vcf", "05.bed");
#get_ref_allele("05.bed", "11.impute", "12.ref.tbl", \@qrys);
#vcf_fill_missing("01.raw.vcf", "12.ref.tbl", "21.vcf");

#call_snp();
#call_indel();
merge_sv_vcfs("51.sv.vcf", "50_sv", \@qrys, $tgt);

sub merge_sv_vcfs {
  my ($fo, $do, $qrys, $tgt) = @_;
  -d $do || make_path($do);
  chdir $do || die "cannot chdir to $do\n";
 
  my $d01 = "01.raw";
  -d $d01  || make_path($d01);
  my @vcfs;
  for my $qry (@qrys) {
    my $di = "$ENV{'misc3'}/$qry\_$tgt/31_sv";
    runCmd("bgzip -c $di/11.sv.vcf > $d01/$qry.vcf.gz");
    runCmd("bcftools index $d01/$qry.vcf.gz");
    push @vcfs, "$d01/$qry.vcf.gz";
  }
  my $vcf_str = join(" ", @vcfs);
  runCmd("bcftools merge -m none -O v $vcf_str -o 02.vcf");

  vcf2tbl("02.vcf", "05.bed");
 
  my $d11 = "11_ref";
  -d $d11 || make_path($d11);
  for my $qry (@$qrys) {
    my $fg = "$ENV{'misc3'}/$qry\_$tgt/23_blat/31.9/gax.bed";
    my $ftb = "$d11/$qry.bp.bed";
    my $ftc = "$d11/$qry.cnt.bed";
    runCmd("bedtools intersect -wao -a 05.bed -b $fg > $ftb");
    runCmd("bedtools intersect -c -a 05.bed -b $fg > $ftc");
    bed_intersect($ftb, $ftc, "$d11/$qry.bed");
    call_ref("$d11/$qry.bed", "$d11/$qry.txt");
  }
  my $fstr = join(" ", map {$d11."/".$_.".txt"} @$qrys);
  runCmd("paste $fstr > 12.ref.tbl");
  vcf_fill_missing("02.vcf", "12.ref.tbl", "15.vcf");

  chdir ".." || die "cannot chdir to ..\n";
  runCmd("bcftools norm -m +any $do/15.vcf -O v -o $fo");
}
sub call_indel {
  runCmd("bcftools view -v indels,mnps,other 21.vcf.gz | bcftools norm -m +any -O v -o 31.indel.vcf");
}
sub call_snp { ### run comp.vcf.R to generate 81.cvg.tbl
  runCmd("bcftools view -v snps -R 81.cvg.tbl -O v 21.vcf.gz | vcf.addref.pl -i - -o 23.snp.vcf");
  runCmd("vcf.stat.py 23.snp.vcf 25.stat.tbl");
}
sub merge_vcfs {
  my ($qrys) = @_;
  -d "00.raw" || make_path("00.raw");
  my @vcfs;
  for my $qry (@qrys) {
    my $di = "$ENV{'misc3'}/$qry\_$tgt/23_blat/31.9";
#    runCmd("bgzip -c $di/vnt.vcf > 00.raw/00.$qry.vcf.gz");
#    runCmd("bcftools index 00.raw/00.$qry.vcf.gz");
    push @vcfs, "00.raw/00.$qry.vcf.gz";
  }
  my $vcf_str = join(" ", @vcfs);
  runCmd("bcftools merge -m none -O u $vcf_str | bcftools norm -m +snps -O v -o 01.raw.vcf");
  runCmd("bgzip -c 01.raw.vcf > 01.raw.vcf.gz");
  runCmd("tabix -p vcf 01.raw.vcf.gz");
}
sub vcf2tbl {
  my ($fi, $fo) = @_;

  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  while( <$fhi> ) {
    chomp;
    next if /(^\#)|(^\s*$)/s;
    my ($chr, $pos, $id, $ref, $alt, $qual, $fil, $info, $fmt, @sams) = 
      split "\t";
    my @alts = split(",", $alt);
    my $reflen = length($ref);
    my @altlens = uniq(map {length($_)} @alts);
    @altlens == 1 || die "$chr:$pos $ref-$alt alts not same lens\n";
    my $altlen = $altlens[0];
    
    if($reflen == 1 && $altlen == 1) {
      print $fho join("\t", $chr, $pos-1, $pos, 'snp')."\n";
    } elsif($reflen > 1 && $altlen == 1) {
      print $fho join("\t", $chr, $pos-1, $pos+$reflen-1, 'del')."\n";
    } elsif($reflen == 1 && $altlen > 1) {
      print $fho join("\t", $chr, $pos-1, $pos+1, 'ins')."\n";
    } else {
      print $fho join("\t", $chr, $pos-1, $pos+$reflen, 'mnp')."\n";
    }
    $alt = $alts[0];
  }
  close $fhi;
  close $fho;
}
sub bed_intersect {
  my ($ftb, $ftc, $fo) = @_;
  open(my $fhb, "<$ftb") or die "cannot read $ftb\n";
  open(my $fhc, "<$ftc") or die "cannot read $ftc\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  while(<$fhc>) {
    chomp;
    my @ps = split "\t";
    my ($chr, $beg, $end, $cnt) = @ps[0,1,2,$#ps];
    if($cnt == 0) {
      my $line = <$fhb>;
      chomp($line);
      my @ps2 = split("\t", $line);
      my ($chr2, $beg2, $end2, $bp2) = @ps2[0,1,2,$#ps2];
      die "$chr:$beg-$end $cnt not 0: $bp2\n" if $bp2 != 0;
      die "sync error: $chr:$beg-$end $cnt\n" if $chr ne $chr2 ||
        $beg != $beg2 || $end != $end2;
      print $fho join("\t", $chr, $beg, $end, @ps[3..$#ps-1], 0, 0)."\n";
      next;
    }

    my $bp = 0;
    for my $i (1..$cnt) {
      my $line = <$fhb>;
      chomp($line);
      my @ps2 = split("\t", $line);
      my ($chr2, $beg2, $end2, $bp2) = @ps2[0,1,2,$#ps2];
      die "sync error: $chr:$beg-$end $cnt\n" if $chr ne $chr2 ||
        $beg != $beg2 || $end != $end2;
      $bp += $bp2;
    }
    print $fho join("\t", $chr, $beg, $end, @ps[3..$#ps-1], $cnt, $bp)."\n";
  }
  close $fhb;
  close $fhc;
  close $fho;
}
sub call_ref {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  while(<$fhi>) {
    chomp;
    my ($chr, $beg, $end, $type, $cnt, $bp) = split "\t";
    my $sta = ($cnt == 1 && $bp == $end - $beg) ? 1 : 0;
    print $fho $sta."\n";
  }
  close $fhi;
  close $fho;
}
sub get_ref_allele {
  my ($fi, $do, $fo, $qrys) = @_;
  -d $do || make_path($do);
  for my $qry (@$qrys) {
    my $fg = "$ENV{'misc3'}/$qry\_$tgt/23_blat/31.9/gax.bed";
    my $ftb = "$do/$qry.bp.bed";
    my $ftc = "$do/$qry.cnt.bed";
    runCmd("bedtools intersect -wao -a 05.bed -b $fg > $ftb");
    runCmd("bedtools intersect -c -a 05.bed -b $fg > $ftc");
    bed_intersect($ftb, $ftc, "$do/$qry.bed");
    call_ref("$do/$qry.bed", "$do/$qry.txt");
  }
  my $fstr = join(" ", map {"$do/".$_.".txt"} @$qrys);
  runCmd("paste $fstr > 12.ref.tbl");
}
sub vcf_fill_missing {
  my ($fi, $fr, $fo) = @_;

  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fhr, "<$fr") or die "cannot read $fr\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  while( <$fhi> ) {
    chomp;
    if(/(^\#)|(^\s*$)/s) {
      print $fho $_."\n";
      next;
    }
    my ($chr, $pos, $id, $ref, $alt, $qual, $fil, $info, $fmt, @gts) = 
      split "\t";
    my $line = <$fhr>;
    chomp($line);
    my @stas = split("\t", $line);
    @gts == @stas || die "fields not equal: $chr:$pos\n";
    for my $i (0..$#gts) {
      my ($gt, $sta) = ($gts[$i], $stas[$i]);
#      die "err: $chr:$pos $gt $sta\n" if $gt eq "1/1" && $sta == 1;
      $gts[$i] = "0/0" if $gt eq "./." && $sta == 1;
    }
    print $fho join("\t", $chr, $pos, $id, $ref, $alt, $qual, $fil, $info, $fmt, @gts)."\n"; 
  }
  close $fhi;
  close $fho;
}

__END__
