package Pileup;
use strict; use Init; use Common; use Path::Class;
use Bio::AlignIO; use Bio::Seq; use Data::Dumper; use Graph; 
use Bio::Factory::FTLocationFactory; use Clone qw/clone/;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use Time::HiRes qw/gettimeofday tv_interval/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT_OK = qw//;
@EXPORT = qw//;

sub hydra2Gff3 {
    my ($inDir, $outDir, $accAry) = rearrange(['indir', 'outdir', 'acc'], @_);
    system("mkdir -p $outDir") unless -d $outDir;
    for my $acc (@$accAry) {
        my $fIn = file($inDir, "$acc.final");
        die "$fIn is not there\n" unless -s $fIn;
        my $fOut = file($outDir, "$acc.gff3");
        bedpe2Gff3($fIn, $fOut, $acc);
    }
}
sub loadHydraGff3 {
    my ($inDir, $accAry, $db) = rearrange(['indir', 'acc', 'db'], @_);
    my $ld = Localdb->new(-db=>$db);
    for my $acc (@$accAry) {
        my $fIn = file($inDir, "$acc.gff3");
        die "$fIn is not there\n" unless -s $fIn;
        $ld->loadGff(-empty=>1, -files=>$fIn);
    }
}


1;
__END__
