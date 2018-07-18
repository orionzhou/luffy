#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.sv.pl - 

=head1 SYNOPSIS
  
  comp.sv.pl [-help] [-qry qry-genome] [-tgt tgt-genome]

  Options:
    -h (--help)   brief help message
    -q (--qry)    qry genome
    -t (--tgt)    tgt genome

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use Tabix;
use Bio::DB::Fasta;
use Common;
use Location;
use Gtb;
use Gal;
use Bed;
use Seq;
use Vcfhead;
use List::Util qw/min max sum/;

my $help_flag;
my ($qry, $tgt) = qw/HM004 HM101/;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "qry|q=s"  => \$qry,
  "tgt|t=s"  => \$tgt,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $dir = "$ENV{'misc3'}/$qry\_$tgt/31_sv";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $tdir = "$ENV{'genome'}/$tgt";
my $qdir = "$ENV{'genome'}/$qry";
my $cdir = "$ENV{'misc3'}/$qry\_$tgt/23_blat";

my $tref = Bio::DB::Fasta->new("$tdir/11_genome.fas");
my $qref = Bio::DB::Fasta->new("$qdir/11_genome.fas");

runCmd("rm *");
cat_var($tdir, $qdir, $cdir, "01.stb");
stb_filt("01.stb", "05.stb");
stb2stx("05.stb", "05.stx");
stb2vcf("05.stb", "10.raw.vcf");
sv_vcf_refine("10.raw.vcf", "11_refine", "11.sv.vcf");
runCmd("vcf.svtbl.py 11.sv.vcf 11.sv.tbl");

##refine_tlc("01.stb", "05.refine.stb");

sub cat_var { # only call SVs in first level alignments
  my ($tdir, $qdir, $dir, $fo) = @_;
  my $fi = "$dir/31.9/idm";
  open(my $fhi, "<$fi") || die "cannot read $fi\n";
  open(my $fho, ">$fo") || die "cannot write $fo\n";
  print $fho join("\t", qw/id tchr tbeg tend tlen tinfo srd
    qchr qbeg qend qlen qinfo/)."\n";

  my $gta = Tabix->new(-data => "$dir/31.5/gax.gz");
  my $gqa = Tabix->new(-data => "$dir/41.5/gax.gz");
  my $gtr = Tabix->new(-data => "$dir/31.9/gax.gz");
  my $gqr = Tabix->new(-data => "$dir/41.9/gax.gz");

  my $gapt = Tabix->new(-data => "$tdir/16.gap.bed.gz");
  my $gapq = Tabix->new(-data => "$qdir/16.gap.bed.gz");

  my $id = 1;
  while( <$fhi> ) {
    chomp;
    next if /(^\#)|(^\s*$)/;
    my ($tid, $tbeg, $tend, $tsrd, $qid, $qbeg, $qend, $qsrd, $cid, $lev) 
      = split "\t";
    $lev == 1 || next;
    $tsrd eq "+" or die "$tid:$tbeg-$tend not +\n";
    my $srd = $qsrd;
    my $tlen = $tend - $tbeg - 1;
    my $qlen = $qend - $qbeg - 1;
   
    my ($flagt, $flagq) = (0, 0);
    if($tlen > 0) {
      my $ary = read_gap($gapt, $tid, $tbeg + 1, $tend - 1);
      $flagt = 1 if @$ary > 0;
    }
    if($qlen > 0) {
      my $ary = read_gap($gapq, $qid, $qbeg + 1, $qend - 1);
      $flagq = 1 if @$ary > 0;
    }
    next if $flagt || $flagq;

    my $tstats = $tlen > 0 ? 
      idm_cat($tid, $tbeg+1, $tend-1, $lev, $gtr, $gta, $gqr) : [];
    my $tinfo = join(",", map {join("-", @$_)} @$tstats);
    my $qstats = $qlen > 0 ? 
      idm_cat($qid, $qbeg+1, $qend-1, $lev, $gqr, $gqa, $gtr) : [];
    my $qinfo = join(",", map {join("-", @$_)} @$qstats);
    print $fho join("\t", $id++, $tid, $tbeg, $tend, $tlen, $tinfo, $srd,
      $qid, $qbeg, $qend, $qlen, $qinfo)."\n";
  }
  close $fhi;
  close $fho;
}
sub idm_cat {
  my ($id, $beg, $end, $lev, $gr, $ga, $gqr) = @_;
  my @stats;
  if($end - $beg + 1 < 30) {
    push @stats, [$beg, $end, 'pav'];
    return \@stats;
  }
  
  my $hr = read_chain($id, $beg, $end, $gr);
  my $lev2 = min( map {$hr->{$_}->[7]} keys(%$hr) );
  my @cids = grep {$hr->{$_}->[7] == $lev2} keys(%$hr);
  my $loct = [];
  for my $cid (@cids) {
    my ($ti, $tb, $te, $qi, $qb, $qe, $len) = @{$hr->{$cid}};
    my $hx = read_chain($qi, $qb, $qe, $gqr);
    my $leno = check_ovlp($ti, $tb, $te, $hx);
    $leno / ($te - $tb + 1) >= 0.8 || next;
    push @stats, [$tb, $te, 'tlc', $qi, $qb, $qe];
    push @$loct, [$tb, $te];
  }
  my $len1 = locAryLen($loct);
  my $len2 = locAryLen(posMerge($loct));
  $len1 == $len2 || die "$id:$beg-$end error\n".Dumper($hr)."\n";
 
  my $locp = [];
  my $hp = read_chain($id, $beg, $end, $ga);
  for my $cid (keys(%$hp)) {
    my ($ti, $tb, $te, $qi, $qb, $qe, $len, $lev2) = @{$hp->{$cid}};
    push @$locp, [$tb, $te];
  }
  $locp = posMerge($locp);
  
#  $beg =~ /^1074/ && die "$id:$beg-$end\n".Dumper(read_gax($ga, $id, $beg, $end)); 
  my ($locd) = posDiff([[$beg, $end]], $locp);
  for (@$locd) {
    my ($tb, $te) = @$_;
    push @stats, [$tb, $te, 'pav'];
  }
  my $locu = posSubtract($locp, $loct);
  for (@$locu) {
    my ($tb, $te) = @$_;
    push @stats, [$tb, $te, 'cnv'];
  }
  return [sort {$a->[0] <=> $b->[0]} @stats];
}
sub stb_filt {
  my ($fi, $fo) = @_;
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  my $t = readTable(-in => $fi, -header => 1);

  my @idxs;
  for my $i (0..$t->lastRow) {
    my ($id, $tchr, $tbeg, $tend, $tlen, $tinfo, $srd,,
      $qchr, $qbeg, $qend, $qlen, $qinfo) = $t->row($i);
    my ($ref, $alt);
    if($tlen > 0) {
      $ref = $tref->seq($tchr, $tbeg, $tend - 1);
    } else {
      $ref = $tref->seq($tchr, $tbeg, $tbeg);
    }
    $alt = substr($ref, 0, 1);
    if($qlen > 0) {
      $alt .= seqret_simple($qref, $qchr, $qbeg+1, $qend-1, $srd); 
    }
    push @idxs, $i if $ref eq $alt;
  }
  $t->delRows(\@idxs);
  print $fho $t->tsv(1);
  close $fho;
}
sub stb2stx {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/id tchr tbeg tend srd qchr qbeg qend type/)."\n";
  while(<$fhi>) {
    chomp;
    /^(id)|(\#)/ && next;
    my ($id, $tchr, $tbeg, $tend, $tlen, $tinfo, $srd, 
      $qchr, $qbeg, $qend, $qlen, $qinfo) = split "\t";
    my @tinfos = split(",", $tinfo);
    my @qinfos = split(",", $qinfo);
    for (@tinfos) {
      my ($tb, $te, $type, $qi, $qb, $qe) = split "-";
      my $line = $type eq "pav" ?
        join("\t", $id, $tchr, $tb, $te, $srd, $qchr, $qbeg, '', 'DEL') :
          $type eq "cnv" ?
        join("\t", $id, $tchr, $tb, $te, $srd, $qchr, $qbeg, '', 'CNL') :
        join("\t", $id, $tchr, $tb, $te, $srd, $qi, $qb, $qe, 'TLC:DEL');
      print $fho $line."\n";
    }
    for (@qinfos) {
      my ($qb, $qe, $type, $ti, $tb, $te) = split "-";
      my $line = $type eq "pav" ?
        join("\t", $id, $tchr, $tbeg, '', $srd, $qchr, $qb, $qe, 'INS') :
          $type eq "cnv" ?
        join("\t", $id, $tchr, $tbeg, '', $srd, $qchr, $qb, $qe, 'CNG') :
        join("\t", $id, $ti, $tb, $te, $srd, $qchr, $qb, $qe, 'TLC:INS');
      print $fho $line."\n";
    }
  }
  close $fhi;
  close $fho;
}
sub refine_tlc {
  my ($fi, $fo) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  stb_tlc2bed($t, "$fo.1.bed");
  runCmd("mergeBed -i $fo.1.bed -c 4 -o collapse > $fo.2.bed");
  
  open(my $fhb, "<$fo.2.bed") or die "cannot read $fo.2.bed\n";
  my @idxs_rm;
  while(<$fhb>) {
    chomp;
    my ($chr, $beg, $end, $idxstr) = split "\t";
    my @idxs = split(",", $idxstr);
    my $h;
    @idxs > 1 || next;
    for my $i (0..$#idxs-1) {
      my ($idx1, $idx2) = ($idxs[$i]-1, $idxs[$i+1]-1);
      my ($flag, $idxk, $idxr) = check_rec_tlc($t, $idx1, $idx2);
      $flag == 1 || next;
      !exists $h->{$idxk} || die "$idxk in >1 rec-tlc\n";
      $h->{$idxk} = $idxr;
      $t->setElm($idxk, "tichr", $t->elm($idxr, 'tichr'));
      $t->setElm($idxk, "tibeg", $t->elm($idxr, 'tibeg'));
      $t->setElm($idxk, "type", 'TLC');
      push @idxs_rm, $idxr;
    }
  }
  close $fhb;

  $t->delRows(\@idxs_rm);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $t->tsv(1);
  close $fho;

  runCmd("rm $fo.*.bed");
}
sub stb_tlc2bed {
  my ($ti, $fo) = @_;
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  for my $i (0..$ti->lastRow) {
    my ($tdchr, $tdbeg, $tdend) = $ti->row($i);
    print $fho join("\t", $tdchr, $tdbeg, $tdend, $i+1)."\n";
  }
  close $fho;
}
sub check_rec_tlc {
  my ($t, $idx1, $idx2) = @_;
  my ($type1, $type2) = map {$t->elm($_, "type")} ($idx1, $idx2);
  my ($idxk, $idxr) = ($type1 eq "qTLC" && $type2 eq "tTLC") ? 
    ($idx2, $idx1) : ($type1 eq "tTLC" && $type2 eq "qTLC") ? 
    ($idx1, $idx2) : (undef, undef);
  return 0 unless defined $idxk;
  my @ary1 = $t->row($idxk);
  my @ary2 = $t->row($idxr);
  my ($ti1, $tb1, $te1, $qi1, $qb1, $qe1) = @ary1[0..2,5..7];
  my ($ti2, $tb2, $te2, $qi2, $qb2, $qe2) = @ary2[0..2,5..7];
  my ($tl1, $ql1) = ($te1 - $tb1 + 1, $qe1 - $qb1 + 1);
  my ($tl2, $ql2) = ($te2 - $tb2 + 1, $qe2 - $qb2 + 1);
  return 0 if $ti1 ne $ti2 || $qi1 ne $qi2;
  my $to = min($te1, $te2) - max($tb1, $tb2) + 1;
  my $qo = min($qe1, $qe2) - max($qb1, $qb2) + 1;
  return 0 if $to/$tl1 < 0.9 || $to/$tl2 < 0.9 ||
    $qo/$ql1 < 0.9 || $qo/$ql2 < 0.9;
  return (1, $idxk, $idxr);
}
sub read_chain {
  my ($chr, $beg, $end, $gax) = @_;
  my $h = {};
  my $ary = read_gax($gax, $chr, $beg, $end);
  for (@$ary) {
    my ($cid, $tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd, $lev) = @$_;
    $h->{$cid} ||= [$tid, $tb, $te, $qid, $qb, $qe, 0, $lev];
    $h->{$cid}->[1] = min($tb, $h->{$cid}->[1]); 
    $h->{$cid}->[2] = max($te, $h->{$cid}->[2]); 
    $h->{$cid}->[4] = min($qb, $h->{$cid}->[4]); 
    $h->{$cid}->[5] = max($qe, $h->{$cid}->[5]); 
    $h->{$cid}->[6] += $te - $tb + 1;
  }
  return $h;
}
sub check_ovlp {
  my ($tid, $tbeg, $tend, $h) = @_;
  my $loc = [];
  for my $cid (keys(%$h)) {
    my ($qi, $qb, $qe, $ti, $tb, $te, $len) = @{$h->{$cid}};
    $ti eq $tid || next;
    push @$loc, [$tb, $te];
  }
  my ($ref, $olen) = posOvlp([[$tbeg, $tend]], $loc);
  return $olen;
}

sub stb2vcf {
  my ($fi, $fo) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $vcfhead."\n";
  print $fho "#".join("\t", @colhead, $qry)."\n";
  for my $i (0..$t->lastRow) {
    my ($id, $tchr, $tbeg, $tend, $tlen, $tinfo, $srd,,
      $qchr, $qbeg, $qend, $qlen, $qinfo) = $t->row($i);
    my ($ref, $alt);
    if($tlen > 0) {
      $ref = $tref->seq($tchr, $tbeg, $tend - 1);
    } else {
      $ref = $tref->seq($tchr, $tbeg, $tbeg);
    }
    $alt = substr($ref, 0, 1);
    if($qlen > 0) {
      $alt .= seqret_simple($qref, $qchr, $qbeg+1, $qend-1, $srd); 
    }
    print $fho join("\t", $tchr, $tbeg, ".", $ref, $alt, 50, '.',
      ".", 'GT', '1/1')."\n";
  }
  close $fho;
}
sub sv_vcf2bed {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  while(<$fhi>) {
    chomp;
    next if /(^\#)|(^\s*$)/s;
    my ($chr, $pos, $id, $ref, $alt, $qual, $fil, $info, $fmt, @sams) = 
      split "\t";
    my @alts = split(",", $alt);
    substr($ref, 0, 1) eq substr($alts[0], 0, 1) || die "chr:$pos $ref $alt\n";
    print $fho join("\t", $chr, $pos - 1, $pos)."\n";
  }
  close $fhi;
  close $fho;
}
sub refine_vcf {
  my ($fi, $fr, $fo) = @_;

  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fhr, "<$fr") or die "cannot read $fr\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  while( <$fhi> ) {
    chomp;
    if(/(^\#)|(^\s*$)/s) {
      print $fho $_."\n";
      next;
    }
    my ($chr, $pos, $id, $ref, $alt, $qual, $fil, $info, $fmt, @gts) = 
      split "\t";
    my $line = <$fhr>;
    chomp($line);
    my @stas = split("\t", $line);
    my ($chr2, $pos2, $snp) = @stas[0,2,6];
    "$chr:$pos" eq "$chr2:$pos2" || die "sync error: $chr:$pos $chr2:$pos2\n";
    
    $alt = $snp.substr($alt, 1) if($snp ne ".");
    print $fho join("\t", $chr, $pos, $id, $ref, $alt, $qual, $fil, $info, $fmt, @gts)."\n"; 
  }
  close $fhi;
  close $fho;
}
sub sv_vcf_refine {
  my ($fi, $do, $fo) = @_;
  -d $do || make_path($do);
  sv_vcf2bed($fi, "$do/01.bed");
  my $f_snp = "$cdir/31.9/snp.bed";
  runCmd("intersectBed -wao -a $do/01.bed -b $f_snp > $do/03.ovlp.bed");
  refine_vcf($fi, "$do/03.ovlp.bed", "$do/05.vcf");
  runCmd("bcftools norm -f $tdir/11_genome.fas -O v -o $fo $do/05.vcf");
}

__END__

