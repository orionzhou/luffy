#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use InitPath;
use Common;
use Bio::SeqIO;
use Gff;
use Seq;
use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
sub goldengate_MergeSnps {
  my $dir = "/export/lab/data/mt35_illumina/variants300/snps/GoldenGate/illumina";
  my @names = qw/Gentz Joelle Marie Monteros random/;
  my $f10 = "$dir/10_snps.tbl";
  my $f_ref = "$ENV{'genome'}/mt_35/41_genome.fa";
  my $f20 = "$dir/20_illumna.csv";

  my $h;
  for my $name (@names) {
    my $fi = "$dir/snps_$name.txt";
    open(FHI, "<$fi");
    my $line = <FHI>;
    
    while(<FHI>) {
      chomp;
      my @ps = split "\t";
      my ($locus, $chr, $pos, $ref, $maf, $cntA, $cntC, $cntG, $cntT, $n) = @ps[0..2,4,5,7..11];
      $chr = "chr$chr" if $chr =~ /^\d+$/;
      
      my $r = {'A'=>$cntA, 'C'=>$cntC, 'G'=>$cntG, 'T'=>$cntT};
      my @alleles = grep {$r->{$_} > 0} keys(%$r);
      die "not biallelic: $chr:$pos\n" unless @alleles == 2;
      my @vars = grep {$_ ne $ref} @alleles;
      die "not biallelic: $name $chr:$pos\n" unless @vars == 1;
      my $var = $vars[0];

      my $stat = [$ref, $var, $maf, $n, $locus, $name];
      $h->{$chr}->{$pos} ||= [];
      push @{$h->{$chr}->{$pos}}, $stat;
    }
    close FHI;
  }

  open(FHO, ">$f20");
  print FHO join("\t", qw/chr pos ref var locus source/)."\n";
  for my $chr (sort(keys(%$h))) {
    my $h2 = $h->{$chr};
    for my $pos (sort {$a<=>$b} keys(%$h2)) {
      my @stats = @{$h2->{$pos}};
      my @refs = map {$_->[0]} @stats;
      die Dumper(@stats) unless uniq(@refs) == 1;
      my @vars = map {$_->[1]} @stats;
      die Dumper(@stats) unless uniq(@vars) == 1;
      my @mafs = map {$_->[2]} @stats;
      die Dumper(@stats) unless uniq(@mafs) == 1;
      my @ns = map {$_->[3]} @stats;
      die Dumper(@stats) unless uniq(@ns) == 1;
      my $locus = join(" ", map {$_->[-2]} @stats);
      my $source = join(" ", map {$_->[-1]} @stats);
      print FHO join("\t", $chr, $pos, $refs[0], $vars[0], $mafs[0], $ns[0], $locus, $source)."\n";
    }
  }
  close FHO;
}
sub goldengate_Convert2Illumina {
  my $dir = "$ENV{'misc3'}/Mt_goldengate";
  my $fi = "$dir/20120718_in.tbl";
  my $fo = "$dir/20120718_illumina.csv";
  my $f_ref = "$ENV{'genome'}/mt_35/41_genome.fa";
  my $sep = ",";
  open(FHO, ">$fo") || die "cannot open $fi for reading\n";
  print FHO join($sep, qw/Locus_Name Sequence Target_Type Genome_Build_Version Chromosome Coordinate Source Source_Version Sequence_Orientation Plus_Minus/)."\n";

  open(FHI, "<$fi");
  my $line = <FHI>;
  while(<FHI>) {
    chomp;
    my @ps = split "\t";
    my ($chr, $pos, $ref, $var, $locus, $src) = @ps;
    $chr = "chr$chr" if $chr =~ /^\d+$/;
    my $loc = [[$pos-60, $pos+60]];
    my $seq = seqRet($loc, $chr, "+", $f_ref);

    my $seq_up = substr($seq, 0, 60);
    my $ref2 = substr($seq, 60, 1);
    my $seq_dw = substr($seq, 61, 60);
    die "not ref: $ref != $ref2\n" unless $ref eq $ref2;
    
    my $seq_str = "$seq_up\[$ref/$var\]$seq_dw";
    print FHO join($sep, $locus, $seq_str, "SNP", "Mt3.5", $chr, $pos, $src, 0, "Forward", "Plus")."\n";
  }
  close FHI;
  close FHO;
}
sub tair9_updates {
  my $dir = "$ENV{'misc2'}/crp.ssp";
  my $fi = "$dir/tair9_updates_raw.tbl";
  my $fo = "$dir/tair9_updates.tbl";
  open(FHO, ">$fo") or die "cannot open $fo for writing\n";
  print FHO join("\t", qw/chr pos1 pos2 nt_old nt_new/)."\n";
  
  my $hb;
  my $t = readTable(-in=>$fi, -header=>1);
  for my $i (1..$t->nofRow) {
    my ($chr, $version, $pos1, $pos2, $type, $nt, $nt_new) = $t->row($i-1);
    if(!defined($hb->{$chr})) {
      print FHO join("\t", $chr, 1, 1, '', '')."\n";
      $hb->{$chr} = 1;
    }

    my $len;
    if($nt =~ /^Nx(\d+)$/) {
      $len = $1;
    } else {
      $len = length($nt);
    }

    if($type eq "substitution") {
      print FHO join("\t", $chr, $pos1, $pos2, $nt, $nt_new)."\n";
    } elsif($type eq "insertion") {
      print FHO join("\t", $chr, $pos1, $pos2, '', '')."\n";
      print FHO join("\t", $chr, $pos1+1, $pos2+1+$len, '', '')."\n";
    } elsif($type eq "deletion") {
      print FHO join("\t", $chr, $pos1, $pos2, '', '')."\n";
      print FHO join("\t", $chr, $pos1+1+$len, $pos2+1, '', '')."\n";
    } else {
      die "unknown type: $type\n";
    }
  }
  close FHO;
}
sub pfam_acc2name {
  my $dir = "$ENV{'data'}/db/pfam";
  my $fp = "$dir/Pfam-A.hmm.info";
  my $fi = "$dir/te.raw.tbl";
  my $fo = "$dir/te.txt";
  my $t = readTable(-in => $fp, -header => 1);
  my %h = map {$t->elm($_, "acc") => $t->elm($_, "id")} (0..$t->lastRow);

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  my $ti = readTable(-in => $fi, -header => 1);
  for my $i (0..$ti->lastRow) {
    my ($id, $name) = $ti->row($i);
    if(! exists $h{$id}) {
      print "not found: $id [$name]\n";
    } else {
      print $fho $h{$id}."\n";
    }
  }
  close $fho;
}
