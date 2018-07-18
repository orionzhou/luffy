#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.liftover.pl - lift over a GTB file using GAX mapping 

=head1 SYNOPSIS
  
  gtb.liftover.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gtb)
    -o (--out)    output file
    -r (--ref)    tgt ref-seq fasta 
    -x (--gax)    aln block file (*.gax.gz)
    -s (--snp)    aln snp file (*.snp.gz)

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
use Bio::DB::Fasta;
use Common;
use Location;
use Seq;
use Gal;
use Tabix;

my ($fi, $fo) = ('') x 2;
my ($fr, $fx, $fs) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "ref|r=s"  => \$fr,
  "gax|x=s"  => \$fx,
  "snp|s=s"  => \$fs,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;
pod2usage(2) if !$fr || !$fx || !$fs;

my ($fhi, $fho);
if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $db = Bio::DB::Fasta->new($fr);
my $gax = Tabix->new(-data => $fx);
my $snp = Tabix->new(-data => $fs);

while(<$fhi>) {
  chomp;
  /^(id)|(\#)/ && next;
  my ($id, $par, $chr, $beg, $end, $srd, 
    $elocs, $ilocs, $clocs, $flocs, $tlocs, $phase, 
    $src, $conf, $cat1, $cat2, $cat3, $note) = split "\t";
  next if $cat1 ne "mRNA";
  die "no CDS-loc for $id\n" unless $clocs;
  my $rcloc = locStr2Ary($clocs); 
  my $riloc = locStr2Ary($ilocs);
  my $cloc = $srd eq "+" ? 
    [ map {[$beg + $_->[0] - 1, $beg + $_->[1] - 1]} @$rcloc ] :
    [ map {[$end - $_->[1] + 1, $end - $_->[0] + 1]} @$rcloc ]; 
  my $iloc = $srd eq "+" ? 
    [ map {[$beg + $_->[0] - 1, $beg + $_->[1] - 1]} @$riloc ] :
    [ map {[$end - $_->[1] + 1, $end - $_->[0] + 1]} @$riloc ];
  my $staloc = $srd eq "+" ?
    [ $beg + $rcloc->[0]->[0] - 1, $beg + $rcloc->[0]->[0] + 1 ] : 
    [ $end - $rcloc->[0]->[0] - 1, $end - $rcloc->[0]->[0] + 1 ];
  my $stoloc = $srd eq "+" ?
    [ $beg + $rcloc->[-1]->[1] - 3, $beg + $rcloc->[-1]->[1] - 1 ] : 
    [ $end - $rcloc->[-1]->[1] + 1, $end - $rcloc->[-1]->[1] + 3 ];
  my $seqsta = seqret_simple($db, $chr, @$staloc, $srd);
  my $seqsto = seqret_simple($db, $chr, @$stoloc, $srd);

  my @locs;
  if($seqsta =~ /^ATG$/i) {
    push @locs, [@$staloc, "start"];
  }
  if($seqsto =~ /^(TAA)|(TGA)|(TAG)$/i) {
    push @locs, [@$stoloc, "stop"];
  }
  for (@$cloc) { push @locs, [@$_, 'cds'] };
  for (@$iloc) { push @locs, [@$_, 'intron'] };
  
  for (@locs) {
    my ($b, $e, $cat) = @$_;
    my $len = $e - $b + 1;
    my $ary = read_gax($gax, $chr, $b, $e);
    for (@$ary) {
      my ($cid, $tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd) = @$_;
      my $len1 = $te - $tb + 1;
      my $snps = read_snp($snp, $tid, $tb, $te);
      my $n_snp = scalar(@$snps);
      my $type;
      if($cat eq "start" || $cat eq "stop") {
        $type = $cat;
        $n_snp == 0 || next;
      } elsif($cat eq "cds") {
        $type = "CDSpart";
      } elsif($cat eq "intron") {
        $type = "intronpart";
      } else {
        die "unknown cat: $cat\n";
      }
      my $nsrd = $tsrd eq $qsrd ? $srd : get_revsrd($srd);
      print $fho join("\t", $qid, 'pz', $type, $qb, $qe, 10, $nsrd, 
        ".", "pri=2;grp=$id;src=T")."\n";
    }
  }
}
close $fhi;


__END__
