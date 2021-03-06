#!/usr/bin/env perl
#
# POD documentation
#-------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  buildhmm.pl - builds HMM profile from one or more protein alignments

=head1 SYNOPSIS
  
  buildhmm.pl [-help] <-aln input alignment-dir> <-hmm output HMM-dir>

  Options:
    -h (--help)   brief help message
    -i (--aln)    directory containing input alignments
    -o (--hmm)    output directory containing built HMMs

=cut
  
#### END of POD documentation.
#-------------------------------------------------------------------------

use strict; 
use FindBin;
use lib "$FindBin::Bin";

use Pod::Usage;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Path qw/make_path remove_tree/;
use Data::Dumper;
use Common; 
use Seq;
use Align;
use List::Util qw/min max sum/;

my ($help_flag, $dir_aln, $dir_hmm);
GetOptions(
  "help|h"   => \$help_flag,
  'aln|i=s'  => \$dir_aln, 
  'hmm|o=s'  => \$dir_hmm, 
) || pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$dir_aln || !$dir_hmm;
pod2usage("cannot open $dir_aln for reading") if ! -d $dir_aln;

sub copy_aln {
  my ($di, $do) = @_;
  if(abs_path($di) eq abs_path($do)) {
    runCmd("mv $di tmp_aln");
    $di = "tmp_aln";
  }
  
  make_path($do) unless -d $do;
  remove_tree($do, {keep_root=>1});

  opendir(DH, $di) or die "cannot open $di: $!\n";
  for my $fname (sort readdir(DH)) {
    next if $fname =~ /^\./;
    my $fi = "$di/$fname";
    my $fo = "$do/$fname";
    open(FHI, "<$fi") or die "cannot open $fi for reading\n";
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    while(<FHI>) {
      chomp;
      my $line = $_;
      if($line =~ /(\/\d+\-\d+)/) {
        my $rep = " " x length($1);
        $line =~ s/\Q$1\E/$rep/;
      }
      print FHO "$line\n";
    }
    close FHO;
  }
  closedir DH;
}
sub get_subgroups {
  my ($dir, $fo) = @_;
  print "Extracting gene family IDs\n";
  
  open(FH, ">$fo") or die "cannot open $fo for writing\n";
  print FH join("\t", qw/family/)."\n";

  opendir(DH, $dir) or die "cannot open $dir: $!\n";
  for my $fname (sort readdir(DH)) {
    next if $fname =~ /^\./;
    my @ids = read_aln_ids( "$dir/$fname" );
    my $fam = $fname;
    $fam =~ s/\.[^.]+$//;
    print FH join("\t", $fam)."\n"; 
  }
  closedir DH;
  close FH;
}

sub trim_aln {
  my ($f_fam, $di, $do) = @_;
  print "trimming alignments\n";
  
  make_path($do) unless -d $do;
  remove_tree($do, {keep_root=>1});
  
  my $f_bin = "trimal";

  my $t = readTable(-in=>$f_fam, -header=>1);
  for my $i (0..$t->nofRow-1) {
    my ($fam) = $t->row($i);
    my $fi = "$di/$fam.aln";
    my $fo = "$do/$fam.aln";
#  my @ids = read_aln_ids($fi);
    runCmd("$f_bin -in $fi -out $fo -gappyout", 0);
#  runCmd("trimal -in $ft1 -out $ft2 -resoverlap 0.75 -seqoverlap 80", 0);
    die "not work for $fam\n" unless -s $fo;
  }
}
sub check_hmm {
  my ($fi) = @_;
  my $f_bin = "hmmemit";
  
  my $tag = 0;
  if( -s $fi && open(JJ, "$f_bin $fi |") ) {
    my @lines;
    while ( <JJ> ){
      chomp;
      push @lines, $_;
    }
    $tag = 1 if $lines[0] =~ /^\>/;
  }
  return $tag;
}
sub aln2hmm {
  my ($f_fam, $di, $do) = @_;
  print "converting alignments to HMMs\n";
  
  make_path($do) unless -d $do;
  remove_tree($do, {keep_root => 1});
  
  my $f_bin = "hmmbuild";
  
  my $t = readTable(-in=>$f_fam, -header=>1);
  for my $i (0..$t->nofRow-1) {
    my ($fam) = $t->row($i);
    my $fi = "$di/$fam.aln";
    my $fo = "$do/$fam.hmm";
    my $ft1 = "$do/$fam.selex";
    aln_fmt_convert($fi, $ft1, 'clustalw', 'selex');
    die "$ft1 is empty\n" unless -s $ft1;
    
    my $ft2 = "$do/$fam.sum";
    while( !check_hmm($fo) ) {
      runCmd("$f_bin --informat selex -o $ft2 $fo $ft1", -1);
    }
    die "$fo is empty\n" unless -s $fo;
    system("rm $ft1 $ft2");
  }
}

sub check_gap {
  my ($fi) = @_;
  my $ai = Bio::AlignIO->new(-file=>"<$fi");
  my $gap = 0;
  while(my $aln = $ai->next_aln()) {
    for my $seq ($aln->each_seq()) {
      if($seq->seq =~ /[\-\.]/) {
        $gap = 1;
        last;
      }
    }
  }
  return $gap;
}
sub get_hmm_stat {
  my ($f_fam, $d_aln, $d_hmm, $fo) = @_;
  print "extracting MSP statistics\n";
  
  my $f_bin1 = "hmmstat";
  my $f_bin2 = "hmmemit";

  my $t = readTable(-in=>$f_fam, -header=>1);
  open(FH, ">$fo");
  print FH join("\t", qw/id nseq length gap consensus/)."\n";
  for my $i (0..$t->nofRow-1) {
    my ($fam) = $t->row($i);
    my $f_hmm = "$d_hmm/$fam.hmm";
    my $f_aln = "$d_aln/$fam.aln";
    
    my $lines = runCmd("$f_bin1 $f_hmm", 2);
    die "cannot get stat for $fam from $f_hmm\n" unless $lines->[-1] =~ /^\s*\d+\s+$fam/;
    my @ps = split " ", $lines->[-1];
    my ($nseq, $len) = @ps[3,5];
    
    $lines = runCmd("$f_bin2 -c $f_hmm", 2);
    die "cannot get con seq for $fam from $f_hmm\n" unless $lines->[0] =~ /^\>$fam/;
    my $seq = join("", @$lines[1..@$lines-1]);

    my $gap = check_gap($f_aln);
    print FH join("\t", $fam, $nseq, $len, $gap, $seq)."\n";
  }
  close FH;
}
sub split_hmm_single {
  my ($fh, $fs, $dir) = @_;
  make_path($dir) unless -d $dir;
  remove_tree($dir, {keep_root => 1});
  
  my $ts = readTable(-in=>$fs, -header=>1);
  for my $i (0..$ts->nofRow-1) {
    my ($id, $nseq, $len) = $ts->row($i);
    my $fo = "$dir/$id.hmm";
    system("hmmfetch $fh $id > $fo");
  }
}

my $dir = $dir_hmm;
my $d11 = "$dir/11_aln";
copy_aln($dir_aln, $d11);
my $f03 = "$dir/03_fam.tbl";
get_subgroups($d11, $f03);
my $d12 = "$dir/12_aln_trim";
trim_aln($f03, $d11, $d12);
my $d15 = "$dir/15_hmm";
aln2hmm($f03, $d12, $d15);
my $f16 = "$dir/16_stat.tbl";
get_hmm_stat($f03, $d12, $d15, $f16);
my $f21 = "$dir/21_all.hmm";
runCmd("cat $d15/*.hmm > $f21");


__END__

