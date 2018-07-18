#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  genomestat.pl - generating genome statistics

=head1 SYNOPSIS
  
  genomestat.pl [-help] [-o organism]

  Options:
    -h (-help)   brief help message
    -o (--org)   organism
    -k (--kmer)  k-mer (default: 60)
    -p (--part)  part number (default: 0)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw/make_path remove_tree/;
use Common;

my ($org, $k) = ('pan4', 60);
my $p = 0;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "org|o=s"  => \$org,
  "kmer|k=i"  => \$k,
  "part|p=i"  => \$p,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$org;

my $dir = "/home/youngn/zhoup/Data/genome/$org/18_stat_k$k";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

#prepare_seq();

#comp_gc();

$p = sprintf "%02d", $p;
#comp_mapp();
#sum_mapp();

#comp_deltag();
#sum_deltag();

sub prepare_seq {
  runCmd("seqtile.pl -i ../11_genome.fas -o 01.fas -step 1 -size $k");
  -d "02_seq" || make_path("02_seq");
  runCmd("split -l 40000000 -d -a 2 01.fas 02_seq/part.");
}

sub comp_mapp {
  -d "15_mapp" || make_path("15_mapp");
  runCmd("bowtie2 -x \$data/db/bowtie2/$org -f -U 02_seq/part.$p \\
    --end-to-end --fast -k 100 -p 16 --reorder | \\
    sammapp.pl -mis 3 | tbl2bed.pl -o 15_mapp/part.$p.bed");
}
sub sum_mapp { 
  runCmd("cat 15_mapp/part.*.bed > 15_mapp.bed");
  runCmd("bedGraphToBigWig 15_mapp.bed ../15.sizes 15_mapp.bw");
}

sub comp_gc {
  runCmd("seqgc -i 01.fas -o 11_gc.tbl");
  runCmd("tbl2bed.pl -i 11_gc.tbl -o 11_gc.bed");
  runCmd("bedGraphToBigWig 11_gc.bed ../15.sizes 11_gc.bw");
}

sub comp_deltag {
  -d "17_deltag" || make_path("17_deltag");
  chdir "17_deltag" || die "cannot chdir to 17_deltag\n";
  runCmd("seqsplit.pl -i ../02_seq/part.$p -n 16 -o part.$p.");
  runCmd("seq 0 15 | xargs -i printf \"%02d\\n\" {} | \\
    parallel hybrid-ss-min -n DNA part.$p.{} -o part.$p.{}");
  runCmd("seq 0 15 | xargs -i printf \"part.$p.%02d\\n\" {} | \\
    parallel -q sh -c \"unafoldo.pl -i {}.dG -s {} | tbl2bed.pl -o {}.bed\"");
  runCmd("seq 0 15 | xargs -i printf \"part.$p.%02d\\n\" {} | \\
    xargs -i rm {}.ct {}.run {}.37.ext {}.37.plot");
  runCmd("cat part.$p.*.bed > part.$p.bed");
  runCmd("rm part.$p.*.*");
  runCmd("rm part.$p.[0-1]*");
}
sub sum_deltag {
  runCmd("cat 17_deltag/part*.bed > 17_deltag.bed");
  runCmd("bedGraphToBigWig 17_deltag.bed ../15.sizes 17_deltag.bw");
}

__END__
