package Vcfhead;
use strict; 
use Data::Dumper;
use List::Util qw/min max sum/; 
use POSIX qw/ceil floor/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/$vcfhead @colhead/;
@EXPORT_OK = qw//;

our $vcfhead = "##fileformat=VCFv4.1
##reference=file:///home/youngn/zhoup/Data/genome/HM101/11_genome.fasta
##contig=<ID=chr1,length=52991155>
##contig=<ID=chr2,length=45729672>
##contig=<ID=chr3,length=55515152>
##contig=<ID=chr4,length=56582383>
##contig=<ID=chr5,length=43630510>
##contig=<ID=chr6,length=35275713>
##contig=<ID=chr7,length=49172423>
##contig=<ID=chr8,length=45569985>
##contig=<ID=chrU,length=29304494>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##FILTER=<ID=q10,Description=\"Quality below 10\">
##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";

our @colhead = qw/CHROM POS ID REF ALT QUAL FILTER INFO FORMAT/;






1;
__END__

