#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use Path::Class;
use Getopt::Long;
use File::Basename;
use File::Path qw/make_path remove_tree/;
use Data::Dumper;

use InitPath;
use Common;
use Bam;
use Medicago;

our $pindel = "$DIR_src/pindel024t";
our $crest = "$DIR_src/crest";
$ENV{PATH} = "$pindel:$crest:$ENV{PATH}";
$ENV{PERL5LIB} = "$crest:$ENV{PERL5LIB}";

my ($opt, $sm) = ('') x 2;
GetOptions('opt|p=s'=>\$opt, 'sm|s=s'=>\$sm);

my $f_ref = "$DIR_genome/Mtruncatula_4.0/11_genome.fa";
my $f_bwa = "$DIR_db/bwa/mt_40";
my $f_rn = "$DIR_misc3/ncgr_fastq/04_fastq_stats.tbl";
my $f_sm = "$DIR_misc3/ncgr_fastq/21_sample.tbl";

my $dir = "$DIR_misc3/hapmap/40_sv";
my $dirI = "$dir/../11_mapping/31_realigned";
my $f_stat = "$dir/../11_mapping/71_stat.tbl";

#run_hydra($dir, $sm);
my $d31 = "$dir/31_pindel";
run_pindel($dirI, $d31, $sm, $f_ref, $f_stat, $f_rn);
my $d33 = "$dir/33_crest";
#run_crest($dirI, $d33, $sm, $f_ref);
#run_genome_strip($dir, $sm);
#run_tigrasv($dir, $sm);

sub run_hydra {
  my ($dir, $sms) = @_;
  my $d01 = "$dir/01_pos_sorted"; 
  my $d02 = "$dir/02_rn_sorted"; 
  my $f03 = "$dir/03_stat.tbl"; 
  my $t = readTable(-in=>$f03, -header=>1);
  
  my $d11 = "$dir/11_stretched";
  my $d13 = "$dir/13_hydra";
  my $d14 = "$dir/14_hydra_filtered";
  
  my $hydra = "\$src/Hydra-Version-0.5.3/bin/hydra";
  for my $sm (@$sms) {
      my $t2 = $t->match_pattern("\$_->[0] eq '$sm'");
      die "cannot find $sm in table\n" unless $t2->nofRow == 1;
      my ($sm, $rns, $is_mld, $is_mno) = map {$t2->elm(0, $_)} qw/sm rns is_mld is_mno/;
      my $rgs_str = join(" ", map {"-g $_"} split(" ", $rns));
      runCmd("bamPickStretched -i $d01/$sm.bam -o $d11/$sm.bed -m $is_mno $rgs_str -r chr5");
      runCmd("$hydra -in $d11/$sm.bed -out $d13/$sm -mld $is_mld -mno $is_mno -ms 5");
  }
}
sub prepare_pindel {
  my ($fi, $fo, $fs, $fr, $sm) = @_;
  open(FHI, "<$fi") or die "cannot read $fi\n";
  while(<FHI>) {
    next if /^id\t/;
    chomp;
    my ($id1, $mate1, $type1, $chr1, $beg1, $end1, $strd1, $cigar1, $mq1, $is1, $rg1, $seq1) = split("\t", $_);
    my $line2 = <FHI>;
    chomp($line2);
    my ($id2, $mate2, $type2, $chr2, $beg2, $end2, $strd2, $cigar2, $mq2, $is2, $rg2, $seq2) = split("\t", $line2);
    die "$id1/$mate1:$chr1 $beg1 $rg1 <> $id2/$mate2:$chr2 $beg2 $rg2\n" if $id1 ne $id2 || $mate1 != 1 || $mate2 != 2;
  } 
  close FHI;
}
sub prepare_pindel_s {
    my ($fi, $fo, $fs, $fr, $sm) = @_;
    
    my $ts = readTable(-in=>$fs, -header=>1);
    my $h = { map {$ts->elm($_, "lb") => $ts->elm($_, "is_median")} (0..$ts->nofRow-1) };
    
    my $tr = readTable(-in=>$fr, -header=>1);
    for my $i (0..$tr->nofRow-1) {
        my ($rn, $lb) = map {$tr->elm($i, $_)} qw/rn lb/;
        next unless exists $h->{$lb};
        $h->{$rn} = $h->{$lb};
    }

    open(FHI, "<$fi") or die "cannot read $fi\n";
    open(FHO, ">$fo") or die "cannot read $fo\n";
    while(<FHI>) {
        next if /^id\t/;
        chomp;
        my ($id1, $mate1, $type1, $chr1, $beg1, $end1, $strd1, $cigar1, $mq1, $is1, $rg1, $seq1) = split("\t", $_);
        my $line2 = <FHI>;
        chomp($line2);
        my ($id2, $mate2, $type2, $chr2, $beg2, $end2, $strd2, $cigar2, $mq2, $is2, $rg2, $seq2) = split("\t", $line2);
        die "$id1/$mate1:$chr1 $beg1 <> $id2/$mate2:$chr2 $beg2\n" if $id1 ne $id2 || $mate1 != 1 || $mate2 != 2;
        
        next if $mq1 eq 0 || $mq2 eq 0;
        my ($id, $type) = ($id1, $type1);

        die "unequal RG: $id1 $rg1 $rg2\n" if $rg1 && $rg2 && $rg1 ne $rg2;
        my $rg = $rg1 ? $rg1 : $rg2;
        die "no is_meadian for $rg\n" unless exists $$h{$rg};
        my $is = $$h{$rg};
        if($type == 1) {
            if( $chr1 ) {
                print FHO "\@$id/$mate2\n";
                print FHO "$seq2\n";
                print FHO join("\t", $strd1, $chr1, $end1, $mq1, $is, $sm)."\n";
            } else {
                print FHO "\@$id/$mate1\n";
                print FHO "$seq1\n";
                print FHO join("\t", $strd2, $chr2, $end2, $mq2, $is, $sm)."\n";
            }
        } elsif($type == 2) {
                print FHO "\@$id/$mate2\n";
                print FHO "$seq2\n";
                print FHO join("\t", $strd1, $chr1, $end1, $mq1, $is, $sm)."\n";
                print FHO "\@$id/$mate1\n";
                print FHO "$seq1\n";
                print FHO join("\t", $strd2, $chr2, $end2, $mq2, $is, $sm)."\n";
        } else {
            die "unknown type: $type at $id\n" unless $type == 3;
            if( $seq2 ) {
                $is = abs($is1);
                print FHO "\@$id/$mate2\n";
                print FHO "$seq2\n";
                print FHO join("\t", $strd1, $chr1, $end1, $mq1, $is, $sm)."\n";
            } else {
                $is = abs($is2);
                print FHO "\@$id/$mate1\n";
                print FHO "$seq1\n";
                print FHO join("\t", $strd2, $chr2, $end2, $mq2, $is, $sm)."\n";
            }
        }
    }
    close FHI;
    close FHO;
}
sub run_pindel {
    my ($dirI, $dir, $sm, $f_ref, $f_stat, $f_rn) = @_;
  
    my $d21 = "$dir/21_orphan";
    make_path($d21) unless -d $d21;
    
#    runCmd("bamPickOrphan -i $dirI/$sm.bam -o $d21/$sm.tbl", 1);
#    runCmd("sort -k1,1 -k2,2n -t \$\'\\t\' -T $DIR_tmp -o $d21/$sm.sorted.tbl $d21/$sm.tbl", 1);
    prepare_pindel("$d21/$sm.sorted.tbl", "$d21/$sm.pindel", $f_stat, $f_rn, $sm);
    
    my $d31 = "$dir/31_predSV";
    make_path($d31) unless -d $d31;
    runCmd("pindel -f $f_ref -p $d21/$sm.pindel -o $d31/$sm -T 4\n", 1);
#    print "pindel2vcf -p $d31/$sm\_D -r $f_ref -R Mt4.0 -d 20130501\n";
#    print "cov_window -i $d11/04_bp.tbl -o $d11/11_cov.tbl -t acc26 -c 1\n";
}
sub run_crest {
    my ($dirI, $dir, $sm, $f_ref) = @_;
    
    my $f_ref_2bit = "$DIR_db/blat/Mtruncatula_4.0.2bit";

    my $gfServer = `cat \$code/pbs/host`;
    $gfServer =~ s/\s*//g;

    my $d11 = "$dir/11_softclip";
    make_path($d11) unless -d $d11;
    runCmd("extractSClip.pl -i $dirI/$sm.bam --ref_genome $f_ref -o $d11", 1);

    my $d31 = "$dir/31_predSV";
    make_path($d31) unless -d $d31;
    runCmd("CREST.pl -f $d11/$sm.bam.cover -d $dirI/$sm.bam \\
        --ref_genome $f_ref --2bitdir $DIR_db/blat -t $f_ref_2bit -o $d31 \\
        --blatserver $gfServer --blatport 1986", 1);
    
    my $f41 = "$dir/41_sum.tbl";
    open(FHO, ">$f41") or die "cannot write to $f41\n";
    print FHO join("\t", qw/acc chr_l pos_l strand_l n_sc_l chr_r pos_r strand_r n_sc_r type cov_l cov_r len_l len_r pct_idty_l pct_reads_nu_l pct_idty_r pct_reads_nu_r beg_con chr_b beg_con_map end_con chr_e end_con_map seq beg_con2 chr_b2 beg_con_map2 end_con2 chr_e2 end_con_map2 seq2/)."\n";
    
    opendir(my $dh, $d31) || die "cannot open $d31\n";
    for my $fn (sort readdir $dh) {
        if($fn =~ /^(HM\d+)\.bam\.predSV/) {
            my $acc = $1;
            my $fi = "$d31/$fn";
            open(FHI, "<$fi") or die "cannot read $fi\n";
            while(<FHI>) {
                chomp;
                my @ps = split("\t", $_);
                if(@ps == 24) {
                    push @ps, ("") x 7;
                } elsif(@ps == 31) {
                } else {
                    die "cannot handle ".@ps." columns\n";
                }
                die join("\t", @ps) if @ps ne 31;
                print FHO join("\t", $acc, @ps)."\n";
            }
            close FHI;
        }
    }
    closedir $dh;
    close FHO;
}
sub run_tigrasv {
    my ($dir, $sms) = @_;
    my $d01 = "$dir/01_pos_sorted"; 
    my $d02 = "$dir/02_rn_sorted"; 
    my $f03 = "$dir/03_stat.tbl"; 
    my $t = readTable(-in=>$f03, -header=>1);
    
    my $d37 = "$dir/37_tigrasv"; 
    print "tigra-sv -R $f_ref $d37/02_sv.tbl -f $d37/01_bamlist.tbl -o $d37/03.fa -r $d37/04_ref.fa\n";
}
sub run_genome_strip {
    my ($dir, $sms) = @_;
    my $d01 = "$dir/01_pos_sorted"; 
    my $d02 = "$dir/02_rn_sorted"; 
    my $f03 = "$dir/03_stat.tbl"; 
    my $t = readTable(-in=>$f03, -header=>1);
    
    my $f_ref_mask = "$DIR_db/bwa/mt_40.mask.fasta";
    my $d35 = "$dir/35_genome_strip"; 
    
    my $chr = "chr5";
    my $lib = "$svtoolkit/lib";
    my $qscript = "$svtoolkit/qscript";
    for my $sm (@$sms) {
        runCmd("java -Xmx3g -cp $lib/gatk/Queue.jar:$lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            org.broadinstitute.sting.queue.QCommandLine \\
            -S $qscript/SVPreprocess.q -S $qscript/SVQScript.q \\
            -gatk $lib/gatk/GenomeAnalysisTK.jar \\
            -cp $lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            -configFile \$conf/svtoolkit/mt.txt -tempDir $DIR_tmp \\
            -md $d35/01_preprocess -R $f_ref -genomeMaskFile $f_ref_mask \\
            -I $d01/$sm.bam -run", 1);
        runCmd("java -Xmx3g -cp $lib/gatk/Queue.jar:$lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            org.broadinstitute.sting.queue.QCommandLine \\
            -S $qscript/SVDiscovery.q -S $qscript/SVQScript.q \\
            -gatk $lib/gatk/GenomeAnalysisTK.jar \\
            -cp $lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            -configFile \$conf/svtoolkit/mt.txt -tempDir $DIR_tmp \\
            -runDirectory $d35 \\
            -md $d35/01_preprocess -R $f_ref -genomeMaskFile $f_ref_mask \\
            -I $d01/$sm.bam -O $d35/02_alt/$sm.vcf \\
            -minimumSize 100 -maximumSize 1000000 \\
            -windowSize 10000000 -windowPadding 10000 -run", 1);
        runCmd("java -Xmx3g -Djava.library.path=$svtoolkit/bwa \\
            -cp $lib/gatk/Queue.jar:$lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            org.broadinstitute.sting.queue.QCommandLine \\
            -S $qscript/SVAltAlign.q -S $qscript/SVQScript.q \\
            -gatk $lib/gatk/GenomeAnalysisTK.jar \\
            -cp $lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            -configFile \$conf/svtoolkit/mt.txt -tempDir $DIR_tmp \\
            -md $d35/01_preprocess -runDirectory $d35 -R $f_ref \\
            -vcf $d35/02_alt/$sm.vcf \\
            -I $d01/$sm.bam -O $d35/03_alt_align/$sm.bam -run", 1);
    }
}




