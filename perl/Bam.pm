package Bam; 
use strict;
use InitPath;
use Common;
use POSIX qw/ceil floor/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/$picard $gatk $svtoolkit
    bam_sort check_bam/;
@EXPORT_OK = qw//;

my $src = $DIR_src;
our $picard = "$src/picard-tools-1.87";
our $gatk = "$src/GenomeAnalysisTK-2.4-7";
our $svtoolkit = "$src/svtoolkit";

sub bam_sort {
  my ($fi, $fo) = @_;
  runCmd("java -Xmx7g -jar $picard/SortSam.jar \\
    VALIDATION_STRINGENCY=LENIENT TMP_DIR=$DIR_tmp \\
    INPUT=$fi OUTPUT=$fo");
}
sub check_bam {
  my ($fi) = @_;
  my $tag = 1;
  $tag = 0 if ! -s $fi;
  my $cmd = "samtools view -H $fi";
  open(JJ, $cmd." 2>&1 |") || die "Failed: $! in \n$cmd\n";
  my ($tag_sq, $tag_pg) = (0, 0);
  while ( <JJ> ){
    chomp;
    $tag = 0 if /EOF marker is absent/;
    $tag = 0 if /fail to open/;
    $tag_sq = 1 if /^\@SQ/;
    $tag_pg = 1 if /^\@PG/;
#    print $_."\n";
  }
  $tag = 0 if $tag_sq+$tag_pg < 2;
#  die join(" ", $tag_sq, $tag_pg, $fi)."\n";
  return $tag;
}
  



1;
__END__
  
