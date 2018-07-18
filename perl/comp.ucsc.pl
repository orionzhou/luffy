#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.ucsc.pl - generate UCSC trackDb.txt

=head1 SYNOPSIS
  
  comp.ucsc.pl [-help] 

  Options:
    -h (--help)   brief help message

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my @qrys = qw/
  HM058 HM125 HM056 HM129 HM060
  HM095 HM185 HM034 HM004 HM050 
  HM023 HM010 HM022 HM324 HM340
/;
my $tgt = "HM101";

my $dir = "$ENV{'misc3'}/comp.ucsc";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my @cols = ('30,144,255', '60,179,113', '255,102,178', '64,64,64');
my @tkeys = qw/track group type shortLabel longLabel bigDataUrl
  visibility priority color/;
my @basetracks = (
  ['gap', 'map', 'bigBed 3', 'Gaps', 'Assembly Gaps', '16.gap.bb', 'squish', 1, '153,0,76'],
  ['gene', 'genes', 'bigBed 12', 'Gene', 'All Genes', '51.gtb.bb', 'pack', 1, '', 'itemRgb=on'],
);

#write_genome_file("genomes.txt");
#build_tracks_ref();
#for (@qrys) { build_tracks_qry($_) };
#create_tgz();

sub write_genome_file {
  my ($fo) = @_;
  
  my $i = 1;
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  for my $org ($tgt, @qrys) {
    my @ps = (
      ['genome', $org],
      ['trackDb', "$org/trackDb.txt"],
      ['groups', 'groups.txt'],
      ['description', $org],
      ['twoBitPath', "$org/11_genome.2bit"],
      ['organism', 'Medicago truncatula'],
      ['defaultPos', 'scf0001:1-10000'],
      ['orderKey', 100 * $i++],
      ['htmlPath', 'mt_description.html']
    );
    if($org eq $tgt) {
      $ps[3]->[1] = "JCVI_Mt4.0";
      $ps[6]->[1] = "chr1:1001-11000";
    }
    for (@ps) {
      print $fho join(" ", @$_)."\n";
    }
    print $fho "\n";
  }
  close $fho;
}
sub build_tracks_ref {
  my $wdir = $tgt;
  -d $wdir || make_path($wdir);
  chdir $wdir || die "cannot chdir to $wdir\n";
  runCmd("rm -rf *");

  my @tracks = @basetracks;
  push @tracks, ['mappability_60mer', 'genomestat', 'bigWig 0 1', 'mapp_60mer', '60mer Mappability', '18_stat_k60/15_mapp.bw', 'full', 1, '64,64,64', 'maxHeightPixels=30:20:10'];
  
  my $dirg = "$ENV{'genome'}/$tgt";
  my $dird = "$ENV{'data'}/db";
  runCmd("cp $dird/blat/$tgt.2bit 11_genome.2bit");
  runCmd("cp $dirg/15.sizes 11_genome.chrom.sizes");
  runCmd("cp $dirg/16.gap.bb .");
  runCmd("cp $dirg/51.gtb.bb .");
  -d "18_stat_k60" || make_path("18_stat_k60");
  runCmd("cp $dirg/18_stat_k60/15_mapp.bw 18_stat_k60/");

  for my $i (0..$#qrys) {
    my $qry = $qrys[$i];
    my $col = $cols[$i%4];
    
    my $diri = "$ENV{'misc3'}/$qry\_$tgt/23_blat/31.9";
    my $diro = "$qry/31.9";
    -d $diro || make_path($diro);
    runCmd("cp $diri/gal.bb $diro/");
    runCmd("cp $diri/gax.bb $diro/");
    runCmd("cp $diri/snp.bb $diro/");
 
    push @tracks, ["$qry\_blat", 'blat', 'bigBed 12', 
      "$qry\_blat", "$qry Blat",
      "$diro/gal.bb", "squish", $i+1, $col];
    push @tracks, ["$qry\_blatx", 'blatx', 'bigBed 6', 
      "$qry\_blatx", "$qry Blatx",
      "$diro/gax.bb", "hide", $i+1, $col];
    push @tracks, ["$qry\_blat_snp", 'blat_snp', 'bigBed 4', 
      "$qry\_blat_snp", "$qry Blat SNP",
      "$diro/snp.bb", "dense", $i+1, $col];
  }

  my $ft = "trackDb.txt";
  open(my $fht, ">$ft") or die "cannot write $ft\n";
  for (@tracks) {
    my @tvals = @$_;
    @tvals >= @tkeys || die Dumper(@tvals)." < ".@tkeys." values\n";
    for my $i (0..$#tkeys) {
      $tvals[$i] eq "" && next;
      print $fht join(" ", $tkeys[$i], $tvals[$i])."\n";
    }
    if(@tvals == @tkeys) {
      print $fht "\n";
      next;
    }
    for my $j (scalar(@tkeys)..$#tvals) {
      my ($tkey, $tval) = split("=", $tvals[$j]);
      print $fht join(" ", $tkey, $tval)."\n";
    }
    print $fht "\n";
  }
  close $fht;
  chdir ".." || die "cannot chdir to ..\n";
}
sub build_tracks_qry {
  my ($qry) = $_;
  my $wdir = $qry; 
  -d $wdir || make_path($wdir);
  chdir $wdir || die "cannot chdir to $wdir\n";
  runCmd("rm -rf *");

  my @tracks = @basetracks;
  
  my $dirg = "$ENV{'genome'}/$qry";
  my $dird = "$ENV{'data'}/db";
  runCmd("cp $dird/blat/$qry.2bit 11_genome.2bit");
  runCmd("cp $dirg/15.sizes 11_genome.chrom.sizes");
  runCmd("cp $dirg/16.gap.bb .");
  runCmd("cp $dirg/51.gtb.bb .");

  my $diri = "$ENV{'misc3'}/$qry\_$tgt/23_blat/41.9";
  my $diro = "$tgt/41.9";
  -d $diro || make_path($diro);
  runCmd("cp $diri/gal.bb $diro/");
  runCmd("cp $diri/gax.bb $diro/");
  runCmd("cp $diri/snp.bb $diro/");
  
  my $col = $cols[0];
  push @tracks, ["$tgt\_blat", 'blat', 'bigBed 12', 
    "$tgt\_blat", "$qry Blat",
    "$diro/gal.bb", "squish", 1, $col];
  push @tracks, ["$tgt\_blatx", 'blatx', 'bigBed 6', 
    "$tgt\_blatx", "$tgt Blatx",
    "$diro/gax.bb", "hide", 1, $col];
  push @tracks, ["$tgt\_blat_snp", 'blat_snp', 'bigBed 4', 
    "$tgt\_blat_snp", "$tgt Blat SNP",
    "$diro/snp.bb", "dense", 1, $col];

  my $ft = "trackDb.txt";
  open(my $fht, ">$ft") or die "cannot write $ft\n";
  for (@tracks) {
    my @tvals = @$_;
    @tvals >= @tkeys || die Dumper(@tvals)." < ".@tkeys." values\n";
    for my $i (0..$#tkeys) {
      $tvals[$i] eq "" && next;
      print $fht join(" ", $tkeys[$i], $tvals[$i])."\n";
    }
    if(@tvals == @tkeys) {
      print $fht "\n";
      next;
    }
    for my $j (scalar(@tkeys)..$#tvals) {
      my ($tkey, $tval) = split("=", $tvals[$j]);
      print $fht join(" ", $tkey, $tval)."\n";
    }
    print $fht "\n";
  }
  close $fht;
  chdir ".." || die "cannot chdir to ..\n";
}
sub create_tgz {
  -s "all.tgz" && runCmd("rm all.tgz");
  runCmd("tar czf all.tgz *");
}
__END__

