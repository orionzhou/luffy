#!/usr/bin/perl -w
#
# POD documentation
#-------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  aln2hmm.pl - convert ALN to HMM

=head1 SYNOPSIS
  
  aln2hmm.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (*.aln) 
    -o (--out)    output file (*.hmm)

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

my ($fi, $fo) = ('') x 2;
my $help_flag;
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my $base = [ split(/\./, $fo) ]->[0];
runCmd("trimAl -in $fi -out $base.1.aln -gappyout");
runCmd("hmmbuild --amino $base.2.hmm $base.1.aln");

sub check_hmm {
  my ($fi) = @_;
  my $f_bin = $ENV{"HMMER"}."/bin/hmmemit";
  die("cannot execute hmmemit: $f_bin not there") unless -s $f_bin;
  
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
  
  my $f_bin = $ENV{"HMMER"}."/bin/hmmbuild";
  die("cannot execute hmmbuild: $f_bin not there") unless -s $f_bin;
  
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
  
  my $f_bin1 = $ENV{"HMMER"}."/bin/hmmstat";
  my $f_bin2 = $ENV{"HMMER"}."/bin/hmmemit";
  die("cannot execute hmmstat: $f_bin1 not there") unless -s $f_bin1;
  die("cannot execute hmmemit: $f_bin2 not there") unless -s $f_bin2;

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


__END__

