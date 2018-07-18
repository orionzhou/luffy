#!/usr/bin/perl -w
use strict; 
use Init; 
use Common; 
use Run; 
use Localdb; 
use Readfile; 
use Path::Class;
use Parser; 
use Gff; 
use Crp; 
use Mapping;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
my ($feDb, $refDb) = ("mt_35", "mt_35");
my $f00 = file($DIR_Misc2, "crp", "06_groups.txt");

my $DIR_Work = dir($DIR_Misc2, "crp", "30_blast");
$DIR_Work = dir($DIR_Misc2, "crp", "31_neighbor");
my $type = "up1k";
my $f01 = file($DIR_Genome, $refDb, "01_seq_$type.fa");
#writeSeqByOpt(-type=>$type, -out=>$f01, -opt=>1, -fedb=>$feDb, -refdb=>$refDb);
my $f02 = file($DIR_Work, "02_ncr_$type.fa");
#writeSeqByOpt(-type=>$type, -out=>$f02, -opt=>2, -fedb=>$feDb, -refdb=>$refDb, -in=>$f00);
my $f03 = file($DIR_Work, "03_blast_out.txt");
#run_blast(-in=>$f02, -out=>$f03, -program=>'blastn', -db=>'mt_35_up1k', -aln=>0);
my $f04 = file($DIR_Work, "04_blast_filtered.txt");
#blastFilter(-in=>$f03, -out=>$f04, -strand=>1);
my $f11 = file($DIR_Work, "11_candidates.txt");
#getCandidates(-in=>$f00, -opt=>1, -fedb=>$feDb, -fblast=>$f04, -out=>$f11);
#getCandidates(-in=>$f00, -opt=>2, -fedb=>$feDb, -out=>$f11);

my $f12 = file($DIR_Work, "12_homology.txt");
my $d13 = dir($DIR_Work, "13_aln");
#findHomology(-in=>$f11, -out=>$f12, -alndir=>$d13, -fedb=>$feDb, -refdb=>$refDb);
my $f14 = file($DIR_Work, "14_homology.txt");
#addExp(-in=>$f12, -out=>$f14, -ds=>'rnaseq', -opt=>[['qId', 2], ['hId', 5]]);

$DIR_Work = dir($DIR_Misc2, "crp", "40_dotplot");
my $f41 = file($DIR_Work, "01_anchors.txt");
my $f42 = file($DIR_Work, "02_loc.txt");
#anchorToLoc(-in=>$f41, -out=>$f42, -fedb=>$feDb, -refdb=>$refDb);
doCombo(-in=>$f42, -dir=>$DIR_Work);


sub mapStatus {
  my ($fi1, $fi2, $fo, $refDb) = rearrange(['in1', 'in2', 'out', 'refdb'], @_);
  my $ld = Localdb->new(-db=>$refDb);
  my $t1 = readTable(-in=>$fi1, -header=>1);
  my $t2 = readTable(-in=>$fi2, -header=>1);
  my $ref1 = [];
  for my $i (0..$t1->nofRow-1) {
    my ($id1, $id2, $fam, $chr, $locS) = map {$t1->elm($i, $_)} qw/id1 id2 family chr location/;
    push @$ref1, [$id1, $id2, $fam, locStr2Obj($locS, $chr)];
  }
  my $ref = {map {$_=>[]} (0..@$ref1-1)};
  for my $i (0..$t2->nofRow-1) {
    my ($fam, $id, $status, $locS) = map {$t2->elm($i, $_)} qw/group hitId status location/;
    my ($chr, $s, $e, $str) = splitLocStr($locS);
    if($chr =~ s/^(ctg)/contig\_/) {
      my $loc1 = Bio::Location::Simple->new(-seq_id=>$chr, -start=>$s, -end=>$e, -strand=>$str);
      my $loc2 = $ld->bac2Chr($loc1);
      ($chr, $s, $e, $str) = map {$loc2->$_} qw/seq_id start end strand/;
    }
    my @idxs = indexes {$_->[-1]->seq_id eq $chr && $_->[-1]->start-10 <= $s && $_->[-1]->end+10 >= $e} @$ref1;
    if(@idxs == 0) {
      print "$id $locS [$status] not identified by new search\n";
    } elsif(@idxs > 1) {
      die "$id $locS > 1 hit: ".join(",", @idxs)."\n";
    } else {
      push @{$ref->{$idxs[0]}}, $i;
    } 
  }
  my $fh = new IO::File $fo, "w";
  print $fh join("\t", qw/id1 id2 family status/)."\n";
  for my $i (0..$t1->nofRow-1) {
    my ($id1, $id2, $fam) = @{$ref1->[$i]}[0..2];
    my $ann = '';
    if(@{$ref->{$i}} > 0) {
      my $idxs2 = $ref->{$i};
      my @anns = map {$t2->elm($_, 'status')} @$idxs2;
      $ann = join(" | ", @anns);
      $ann ||= "blank";
    }
    print $fh join("\t", $id1, $id2, $fam, $ann)."\n";
  }
}
sub correct4RedGenes {
  my ($fi1, $fi2, $fo) = rearrange(['in1', 'in2', 'out'], @_);
  my $t1 = readTable(-in=>$fi1, -header=>1);
  my $t2 = readTable(-in=>$fi2, -header=>0);
  my $gh = { map {$_=>1} $t2->col("col1") };
  my $fh = new IO::File $fo, "w";
  print $fh join("\t", qw/id1 id2 is_redundant/)."\n";
  for my $i (0..$t1->nofRow-1) {
    my $id = $t1->elm($i, "hit_id");
    $id =~ s/\.1$// if $id =~ /^contig/;
    my $red = '';
    $red = 1 if exists $gh->{$id};
    print $fh join("\t", $t1->elm($i, "id1"), $t1->elm($i, "id2"), $red)."\n";
  }
}

sub alignCrp {
  my ($fi, $f_seq, $f_hmm, $dirO) = rearrange(['fi', 'fseq', 'fhmm', 'out'], @_);
  system("mkdir -p $dirO") unless -d $dirO;
  my $t = readTable(-in=>$fi, -header=>1);
  my $h = readSeq($f_seq, 2);
  my $h_hmm = readSeq($f_hmm, 2);

  my $h2;
  for my $i (0..$t->nofRow-1) {
    my ($id, $fam) = map {$t->elm($i, $_)} qw/parent cat3/;
    die "$id: no seq\n" unless exists $h->{$id};
    my $seq = $h->{$id};
    $h2->{$fam} ||= [];
    push @{$h2->{$fam}}, $seq;
  }

  for my $fam (keys %$h2) {
    my $seqs = $h2->{$fam};
    die "$fam: no consensus seq\n" unless exists $h_hmm->{"$fam.trim-consensus"};
    my $seq_con = $h_hmm->{"$fam.trim-consensus"};
    run_clustalw(-seqs=>[$seq_con, @$seqs], -out=>"$dirO/$fam.aln", -type=>"protein");
#    run_tcoffee(-seqs=>[$seq_con, @$seqs], -out=>"$dirO/$fam.aln");
  }
}

sub get_core_seq {
  my ($fi, $fo, $f_seq) = @_;
  my $t = readTable(-in=>$fi, -header=>1);
  my (@seqs_cds, @seqs_pro);
  for my $i (0..$t->nofRow-1) {
    my ($id, $chr, $locStr) = map {$t->elm($i, $_)} qw/parent chr note/;
    my $loc = locStr2Obj($locStr, $chr);
    my $seqStr = seqRet($loc, $f_seq);
    my $seq_cds = Bio::Seq->new(-id=>$id, -seq=>$seqStr);
    my $seq_pro = $seq_cds->translate();
    push @seqs_pro, $seq_pro;
  }
  writeSeq(\@seqs_pro, $fo);
}

sub findHomology {
  my ($fIn, $fOut, $alnDir, $feDb, $refDb, $param) = 
    rearrange(['in', 'out', 'alndir', 'fedb', 'refdb', 'param'], @_);
  system("mkdir -p $alnDir") unless -d $alnDir;
  my $t = readTable(-in=>$fIn, -header=>1);
  my $fOutH = new IO::File $fOut, "w";
  my ($qryLd, $hitLd) = map {Localdb->new(-db=>$_)} ($feDb, $refDb);
  my ($cntHit, $cnt) = (1, 0);
  print $fOutH join("\t", qw/qGroup qId dist hId hGroup 
    up200+in100[pct_idty] upLen_qry upLen_hit aln_id/)."\n";
  for my $i (0..$t->nofRow-1 ) {
    my ($qGroup, $qId, $dist, $hId, $hGroup) = $t->row($i);
    my ($seqQ, $lenQ1, $lenQ2) = $qryLd->retLeadSeq(-id=>$qId, -refdb=>$refDb, -opt=>2);
    my ($seqH, $lenH1, $lenH2) = $hitLd->retLeadSeq(-id=>$hId, -refdb=>$refDb, -opt=>2);
    my ($tag, $pre) = (1, sprintf("%03d", $cntHit));
    last if ++$cnt < 1;
    my $rst;
    $tag = 0 if $lenQ1 == 0 || $lenH1 == 0;
    my $f1 = file($alnDir, "$pre.needle");
    my $ref = align(-seqs=>[$seqQ, $seqH], -out=>$f1, -type=>"DNA", -program=>'needle');
    my ($startQ, $endQ) = map {$ref->{$qId}->{$_}} qw/start end/;
    my ($startH, $endH) = map {$ref->{$hId}->{$_}} qw/start end/;
    die Dumper($ref) if !$startQ || !$startH;
    my $aln_len = $ref->{$qId}->{length};
    my $pct_idty = sprintf("%.03f", $ref->{matches} / $ref->{aln_len});
    $tag = 0 if $startQ > $lenQ1 || $startH > $lenH1;
    $tag = 0 if $endQ < $lenQ1 || $endH < $lenH1;
    $tag = 0 if $aln_len < 100;
    $tag = 0 if $pct_idty < 0.5;
    if($tag == 1) {
      my @tmp = ($pct_idty);
      my $f2 = file($alnDir, "$pre.png");
      ($seqQ, $lenQ1, $lenQ2) = $qryLd->retLeadSeq(-id=>$qId, -refdb=>$refDb, -opt=>1);
      ($seqH, $lenH1, $lenH2) = $hitLd->retLeadSeq(-id=>$hId, -refdb=>$refDb, -opt=>1);
      push @tmp, ($lenQ1, $lenH1);
      run_dotplot(-seqs=>[$seqQ, $seqH], -out=>$f2);
      my $f3 = file($alnDir, "$pre.water");
      $ref = align(-seqs=>[$seqQ, $seqH], -out=>$f3, -type=>"DNA", -program=>'water');
      print $fOutH join("\t", $qGroup, $qId, $dist, $hId, $hGroup, @tmp, $pre)."\n";
      my @types = qw/CDS protein/;
      $cntHit ++;
    } else {
      system("rm ".file($alnDir, "$pre*"));
    }
  }
}
