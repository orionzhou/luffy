#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.syn.pl - extract synonymous positions from a Gtb file

=head1 SYNOPSIS
  
  gtb.syn.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gtb)
    -o (--out)    output file (Tbl)

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
use Data::Dumper;
use Common;
use Location;
use Seq;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my ($fhi, $fho);
if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read file $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $fcod = "$ENV{'misc2'}/codon.degeneracy/codon.degeneracy.tbl";
my $tcod = readTable(-in => $fcod, -header => 1);

my $hc;
for my $i (0..$tcod->lastRow) {
  my @ps = $tcod->row($i);
  my ($codon, $d1, $d2, $d3) = @ps[0,7..9];
  $hc->{$codon} = [$d1, $d2, $d3];
}
#print Dumper($hc)."\n";

my $seqdb = "$ENV{'genome'}/HM101/11_genome.fas";
#print $fho join("\t", qw/id seq degen/)."\n";

while( <$fhi> ) {
  chomp;
  /(^id)|(^\#)|(^\s*$)/i && next;
  my $ps = [ split("\t", $_, -1) ];
  @$ps >= 18 || die "not 19 fileds:\n$_\n";
  my ($id, $par, $chr, $beg, $end, $srd, 
    $locES, $locIS, $locCS, $loc5S, $loc3S, $phaseS,
    $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
  $cat1 eq "mRNA" || next;
  die "no CDS-loc for $id\n" unless $locCS;
  $locCS || die "no CDS for $id\n";
  my $rloc = locStr2Ary($locCS);
  my $loc = $srd eq "-" ? [map {[$end-$_->[1]+1, $end-$_->[0]+1]} @$rloc] : 
    [map {[$beg+$_->[0]-1, $beg+$_->[1]-1]} @$rloc];
  my $len = locAryLen($loc);
  my @poss = loc2pos($loc);

  my $seq = seqRet($loc, $chr, $srd, $seqdb);
  $len == length($seq) || die "unequal length: $id $len <> ".length($seq)."\n";
  my @deg = get_codon_degen($seq, $hc);
  @deg = reverse(@deg) if $srd eq "-";
  for my $i (0..$#poss) {
    print $fho join("\t", $chr, $poss[$i], $deg[$i])."\n";
  }
}
close $fhi;
close $fho;

sub loc2pos {
  my ($loc) = @_;
  my @loce;
  for (sort {$a->[0] <=> $b->[0]} @$loc) {
    my ($a, $b) = @$_;
    push @loce, ($a..$b);
  }
  return @loce;
}
sub get_codon_degen {
  my ($seq, $hc) = @_;
  my @deg;
  my $len = length($seq);
  my $ncodon = int($len / 3);
  for my $i (1..$ncodon) {
    my $codon = substr($seq, ($i-1)*3, 3);
    if(exists($hc->{$codon})) {
      push @deg, @{$hc->{$codon}};
    } else {
      push @deg, ("-") x 3;
    }
  }
  if($len % 3 > 0) {
    push @deg, ("-") x ($len % 3);
  }
  @deg == $len || die "codon degen err\n";
  return @deg;
}

__END__
