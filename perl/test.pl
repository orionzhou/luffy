#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use InitPath;
use Common;
use Gff;
use Gtb;
use Seq;
use Align;
use Data::Dumper;
use Time::HiRes qw/gettimeofday tv_interval/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
sub test_time {
  my $t0 = [gettimeofday];
  my $t1 = [gettimeofday];
  printf "***%.01f min\n", tv_interval($t0, $t1) / 60;
}
sub gmSnp {
    my $dir = "$DIR_data/misc2/gm_snp2";
    my $paramR = [
        {fname=>'01_LG_G_1_10.txt', head=>1, left=>1},
        {fname=>'01_LG_G_1_10_noSSRs.txt', head=>3, left=>1},
        {fname=>'summer07_F.txt', head=>3, left=>1},
        {fname=>'summer07_G.txt', head=>3, left=>1},
        {fname=>'summer07_A2.txt', head=>3, left=>1},
        {fname=>'summer07_J.txt', head=>3, left=>1},
        {fname=>'summer07_O.txt', head=>3, left=>1},
        {fname=>'01_LG_G_oct08.txt', head=>2, left=>1} ];
    for my $param (@$paramR[2]) {
        my $fi = $dir."/".$param->{fname};
        #writeGmSnp($fi, $param);
    }
}
sub testDeepMerge {
    my $a = [ [[6,15]], [[25,35]], [[1,10], [21,30]], [[50,60], [70,80]], [[60,64]], [[90,100]] ];
    my $b = posMergeDeep($a);
    for (@$b) {
        my ($locA, $idxs) = @$_;
        my $locStr = join(" ", map {join("_", @$_)} @$locA);
        print join("\t", $locStr, join(",", @$idxs))."\n";
    }
}
sub test_pg_connect {
    use DBI;
    my $pg_db = 'chado';
    my ($pg_host, $pg_user, $pg_pw) = ($ENV{'PGSQLH'}, $ENV{'PGSQLU'}, $ENV{'PGSQLP'});
    my $dbh = DBI->connect("dbi:Pg:dbname=$pg_db;host=$pg_host", $pg_user, $pg_pw)
        or die "cannot connect to pgsql\n";
}
sub convId2Url {
    my $dir = "$DIR_in/seminar_March_7";
    my @fns = qw/frog_a bovine_a human_a mouse_a frog_b bovine_b human_b mouse_b/;
    my $seqNH = Bio::SeqIO->new(-file=>">$dir/hemoglobin_na.fa", -format=>"fasta");
    my $seqAH = Bio::SeqIO->new(-file=>">$dir/hemoglobin_aa.fa", -format=>"fasta");
    for my $fn (@fns) {
        my $f_gb = file($dir, "01_seqs/$fn.gb");
        my $seqH = Bio::SeqIO->new(-file=>$f_gb, -format=>"genbank");
        my $seq = $seqH->next_seq();
        print join("\t", $seq->id, $seq->length)."\n";
        my @fes_cds = grep {$_->primary_tag eq "CDS"} $seq->get_SeqFeatures();
        die " not 1 CDS in $fn\n" unless @fes_cds == 1;
        my $fe = $fes_cds[0];
        my $seqStr = $seq->subseq($fe->start, $fe->end);
        my $seqCds = Bio::Seq->new(-id=>$fn, -seq=>$seqStr);
        $seqCds = $seqCds->revcom() if $fe->strand == -1;
        $seqNH->write_seq($seqCds);
        print "\t".join("\t", $seqCds->id, $seqCds->length)."\n";
        my $seqPro = $seqCds->translate();
        $seqAH->write_seq($seqPro);
        print "\t".join("\t", $seqPro->id, $seqPro->length)."\n";
    }
}
sub cp_fq_scratch {
    my $fi = "$DIR_misc3/hapmap/09_fastq.tbl";
    my $do = "/project/scratch/zhoup/data_for_taehyun/fastq";
    my $t = readTable(-in=>$fi, -header=>1);
    my %h_sm = map { $_=>1 } qw/HM004 HM010 HM050 HM056 HM101/;
    for my $i (0..$t->nofRow-1) {
        my ($idx, $sm, $rg, $lb, $id_lb, $id_run, $pl, $rl, $pi, $f1a, $f1b) = $t->row($i);
        next unless exists $h_sm{$sm};
        my $fo1 = "$do/$rg.1.fq.gz";
        my $fo2 = "$do/$rg.2.fq.gz";
        runCmd("cp $f1a $fo1");
        runCmd("cp $f1b $fo2");
    }
}
sub fasta_mask {
#my $fi = "/project/youngn/zhoup/Data/misc3/spada/Zmays/01_genome/01_refseq.fa";
#my $fo = "/project/youngn/zhoup/Data/misc3/spada/Zmays_masked/01_genome/01_refseq.fa";
    my ($fi, $fo) = @_;
    open(FHI, "<$fi");
    open(FHO, ">$fo");
    while(<FHI>) {
        chomp;
        if(/^\>/) {
            print FHO $_."\n";
            print $_."\n";
        } else {
            $_ =~ s/[atcgn]/N/g;
            print FHO $_."\n";
        }
    }
    close FHI;
    close FHO;
}
sub double_genome {
#double_genome($fi, $fo);
    my ($fi, $fo) = @_;
    my $seqHI = Bio::SeqIO->new(-file=>"<$fi", -format=>'fasta');
    my $seqHO = Bio::SeqIO->new(-file=>">$fo", -format=>"fasta");
    while(my $seq = $seqHI->next_seq()) {
        print join("\t", $seq->id, $seq->length)."\n";
        $seqHO->write_seq($seq);
        my $seqN = Bio::Seq->new(-id=>$seq->id."_N", -seq=>"N"x$seq->length);
        $seqHO->write_seq($seqN);
    }
    $seqHI->close();
    $seqHO->close();
}
sub vcf_parse {
  use Vcf;
  my $dir = "$DIR_misc3/hapmap_mt40/40_sv";
  my $fi = "$dir/HM056_chr5_SI.vcf";
  my $fo = "$dir/HM056_chr5_SI.tbl";
  open(FHO, ">$fo") or die "cannot write $fo\n";
  my $vcf = Vcf->new(file=>$fi);
  $vcf->parse_header();

  my $header_printed=0;

  while (my $x=$vcf->next_data_hash()) {
    if ( !$header_printed ) {
        print FHO "chr\tpos\ttype\tend\tsvlen\tntlen\thomlen\thomseq\tref\talt";
        for my $col (sort keys %{$$x{gtypes}}) {
            print FHO "\t$col";
        }
        print FHO "\n";

        $header_printed = 1;
    }

    my $y = $$x{INFO};
    my ($ntlen, $homseq) = (0, "");
    $ntlen = $$y{NTLEN} if exists $$y{NTLEN};
    $homseq = $$y{HOMSEQ} if exists $$y{HOMSEQ};
    print FHO join("\t", $$x{CHROM}, $$x{POS}, $$y{SVTYPE}, $$y{END}, $$y{SVLEN}, 
        $ntlen, $$y{HOMLEN}, $homseq, $$x{REF}, $$x{ALT}->[0]);
    for my $col (sort keys %{$$x{gtypes}}) {
        my ($al1,$sep,$al2) = exists($$x{gtypes}{$col}{GT}) ? $vcf->parse_alleles($x,$col) : ('.','/','.');
        my $gt = $al1.'/'.$al2;
        print FHO "\t".$gt;
    }
    print FHO "\n";
  }
}
my $dir = "/home/youngn/zhoup/Data/misc3/HM034_HM101/23_blat";
sv2bed("$dir/31.9/sv.ins.tbl", "$dir/lins.bed");
sv2bed("$dir/31.9/sv.del.tbl", "$dir/ldel.bed");
sub sv2bed {
  my ($fi, $fo) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  open(my $fho, ">$fo") or die "cannot write $fo\n";

  for my $i (0..$t->lastRow) {
    my ($chr, $beg, $end, $loc) = $t->row($i);
    my ($id, $b, $e) = parse_locstr($loc);
    my $len = $e - $b + 1;
    $len >= 200 || next;
    print $fho join("\t", $chr, $beg-1, $end)."\n";
  }
  close $fho;
}
sub tlc2tbl {
  my ($fi, $fo) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  
  my $j = 1;
  for my $i (0..$t->lastRow) {
    my ($ti, $tb, $te, $qi, $qb, $qe, $type, $tiloc, $tdloc) = $t->row($i);
    $type eq "r" || next;
    my ($chr1, $beg1, $end1) = parse_locstr($tiloc);
    my ($chr2, $beg2, $end2) = parse_locstr($tdloc);
    my $len = $te - $tb + 1;
    $len >= 100 || next;
    $beg1 -= int($len / 2);
    $beg2 -= int($len / 2);
    $end1 += int($len / 2);
    $end2 += int($len / 2);
    print $fho join("\t", $j, $chr1, $beg1, $end1)."\n";
    print $fho join("\t", $j, $chr2, $beg2, $end2)."\n";
  }
  close $fho;
}

