#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use File::Path qw/make_path remove_tree/;
use Eutils;
use Common;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $org = "Athaliana";
$org = "Mtruncatula_3.5";
my $dir = "/home/youngn/zhoup/Data/misc2/spada.rnaseq/$org";

my $f01 = "$dir/01_acc.txt";
my $f02 = "$dir/02_uid.tbl";
#get_sra_uid($f01, $f02);
my $f03 = "$dir/03_info.tbl";
#get_sra_info($f02, $f03);

my $d05 = "$dir/05_models";
#gtb_compare.pl -q 03_crp_gs.gtb -t 01_gene.gtb -c 05_gene_gs.tbl -o 05_gene_gs.gtb
#gtb_compare.pl -q 11_crp_spada.gtb -t 05_gene_gs.gtb -c 15.tbl -o 15.gtb
#gtb_conv.pl -i 15.gtb -o 15.gff

my $f_gff = "$d05/15.gff";
my $d11 = "$dir/11_reads";
#download_sra($f03, $d11);

my $d21 = "$dir/21_tophat";
my $f_bt2 = "$DIR_db/bowtie/$pre";
#run_tophat($f03, $d11, $d21, $f_bt2, $f_gff);

my $d23 = "$dir/23_cufflinks";
my $f23 = "$d23/isoforms.fpkm_tracking";
run_cufflinks($d21, $f_gff, $d23);

my $d31 = "$dir/31_rnaseq";
make_path($d31) unless -d $d31;
my $f05_05 = "$d05/05_gene_gs.tbl";
my $f05_15 = "$d05/15.tbl";
get_exp_support($f23, $f05_05, "$d31/05_gs.tbl");
get_exp_support($f23, $f05_15, "$d31/15_spada.tbl");

sub get_exp_support {
  my ($ff, $fi, $fo) = @_;
  
  my $tf = readTable(-in=>$ff, -header=>1);
  my $hf;
  for my $i (0..$tf->nofRow-1) {
    my ($id, $class, $idN, $idG, $name, $idT, $locus, $len, $cov, $fpkm, $fpkm_lo, $fpkm_hi, $fpkm_status) = $tf->row($i);
    $hf->{$id} = [$cov, $fpkm, $fpkm_status];
  }

  my $ti = readTable(-in=>$fi, -header=>1);
  open(FHO, ">$fo") or die "cannot write to $fo\n";
  print FHO join("\t", qw/id cov fpkm fpkm_status/)."\n";
  for my $i (0..$ti->nofRow-1) {
    my ($idQ, $tag, $idT) = $ti->row($i);
    my $stat = ["", "", ""];
    if($tag == 1) {
      if(exists $hf->{$idT}) {
        $stat = $hf->{$idT};
      } else {
        print "no exp for $idT\n";
      }
    } elsif($tag == 10) {
      next;
    } else {
      if(exists $hf->{$idQ}) {
        $stat = $hf->{$idQ};
      } else {
        print "no exp for $idQ\n";
      }
    }
    print FHO join("\t", $idQ, @$stat)."\n";
  }
  close FHO;
}

sub get_sra_uid {
  my ($fi, $fo) = @_;
  my @accs;
  open(FHI, "<$fi") or die "cannot read $fi\n";
  while(<FHI>) {
    chomp;
    next if /(^\#)|(^\s*$)/;
    push @accs, $_;
  }
  close FHI;
  
  open(FHO, ">$fo") or die "cannot write $fo\n";
  print FHO join("\t", qw/acc uid/)."\n";
  for my $acc (@accs) {
    my $fac = Bio::DB::EUtilities->new(-eutil=>'esearch', -email=>'zhoupenggeni@gmail.com', -db=>'sra', -term=>$acc);
    my @ids = $fac->get_ids;
    if($fac->get_count == 0) {
      die "$acc: not found\n";
    } elsif($fac->get_count > 1) {
      die "$acc: >1 found\n";
    } else {
      print FHO join("\t", $acc, $ids[0])."\n";
    }
  }
  close FHO;
}
sub get_sra_info {
  my ($fi, $fo) = @_;
  my $t = readTable(-in=>$fi, -header=>1);
  my @uids = uniq($t->col("uid"));
  
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/uid n_run runs title platform experiment 
    study taxid lib lib_strategy lib_layout/)."\n";
  my $fac = Bio::DB::EUtilities->new(
    -eutil=>'esummary', 
    -email=>'zhoupenggeni@gmail.com',
    -db=>'sra', -id=>\@uids);
  my @its;
  my $xml;
  my $pr = XML::Parser->new( Style => 'Debug' );
  while(my $ds = $fac->next_DocSum) {
    my $uid = $ds->get_id;
    @its = $ds->get_Items_by_name("ExpXml");
    $xml = $its[0]->get_content;
    my ($title, $platform, $n_run, $exp, $study, $taxid) = ("") x 6;
    my ($lib, $lib_strategy, $lib_layout) = ("") x 3;
    $title = $1 if $xml =~ /<Title>([^<]+)<\/Title>/;
    $platform = $1 if $xml =~ /<Platform instrument_model="([\w\s]+)"/;
    $n_run = $1 if $xml =~ /<Statistics total_runs="([\w\s]+)"/;
    $exp = $1 if $xml =~ /<Experiment[^>]+name="([^>]+)"/;
    $study = $1 if $xml =~ /<Study[^>]+name="([^>]+)"/;
    $taxid = $1 if $xml =~ /<Organism taxid="(\d+)"/;
    $lib = $1 if $xml =~ /<LIBRARY_NAME[^>]*>([^<]+)<\/LIBRARY_NAME>/;
    $lib_strategy = $1 if $xml =~ /<LIBRARY_STRATEGY[^>]*>([^<]+)<\/LIBRARY_STRATEGY>/;
    $lib_layout = $1 if $xml =~ /<LIBRARY_LAYOUT[^<]+<(\w+)/;

    @its = $ds->get_Items_by_name("Runs");
    $xml = $its[0]->get_content;
    my @runs;
    while($xml =~ /Run acc="(\w+)"[^>]+is_public="true"/ig) { push @runs, $1; }
    my $runs = join(" ", @runs);
    die "not $n_run runs: $runs\n" unless $n_run == @runs;

    print FHO join("\t", $uid, $n_run, $runs, $title, $platform,
      $exp, $study, $taxid, $lib, $lib_strategy, $lib_layout)."\n";
  }
  close FHO;
}
sub download_sra {
  my ($fi, $do) = @_;
  my $t = readTable(-in=>$fi, -header=>1);
  for my $i (0..$t->nofRow-1) {
#      next if $i < 29;
    my ($n_run, $runs, $lib_layout) = map {$t->elm($i, $_)} qw/n_run runs lib_layout/;
    my @runs = split(" ", $runs);
    die "not $n_run runs: $runs\n" unless $n_run == @runs;
    my $tag = $lib_layout =~ /PAIRED/i ? " --split-3" : "";
    for my $run (@runs) {
      runCmd("fastq-dump --origfmt --gzip --outdir $do $run", 1);
    }
  }
}
sub run_tophat {
  my ($fi, $di, $do, $f_bt2, $f_gff) = @_;
  my $t = readTable(-in=>$fi, -header=>1);
  my @fqs;
  for my $i (0..$t->nofRow-1) {
    my ($n_run, $runs, $lib_layout) = map {$t->elm($i, $_)} qw/n_run runs lib_layout/;
    my @runs = split(" ", $runs);
    for my $run (@runs) {
      my $fq = "$di/$run.fastq.gz";
      die "$run is not there:$fq\n" unless -s $fq;
      push @fqs, $fq;
    }
  }
  printf "%d runs\n", scalar(@fqs);
  my $fi_str = join(",", @fqs);
  my $cmd = "tophat -p 16 --max-intron-length 10000 -o $do -G $f_gff $f_bt2 $fi_str";
  runCmd($cmd, 1);
  runCmd("samtools index $do/accepted_hits.bam", 1);
}
sub run_cufflinks {
  my ($di, $f_gff, $do) = @_;
  my $fi = "$di/accepted_hits.bam";
  make_path($do) unless -d $do;
  runCmd("cufflinks -p 16 -g $f_gff $fi -o $do", 1);
}


