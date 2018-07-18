package Gatk; 
use strict;
use Common;
use POSIX qw/ceil floor/;
use List::Util qw/min max sum/;
use File::Path qw/make_path remove_tree/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/$picard $gatk/;
@EXPORT_OK = qw//;

our $picard = "$ENV{'src'}/picard-tools-1.121";
our $gatk = "$ENV{'src'}/GenomeAnalysisTK.jar";

my $gatk_temp = "/scratch/zhoup/gatk_temp";
-d $gatk_temp || make_path($gatk_temp);

$ENV{"_JAVA_OPTIONS"} = "-Djava.io.tmpdir=$gatk_temp";


1;
__END__
  
