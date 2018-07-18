#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 SYNOPSIS
  
  gtb.merge.pl [-help] [-a input-A] [-b input-B] [-out output-file]

  Options:
    -h (--help)   brief help message
    -a (--ina)    input file A (Gtb, reduce)
    -b (--inb)    input file B (Gtb, keep)
    -o (--out)    output file (Gtb)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Location;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fa, $fb, $fo) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "ina|a=s"  => \$fa,
  "inb|b=s"  => \$fb,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fa || !$fb || !$fo;

gtb2bed($fa, "$fo.a.bed", 2);
gtb2bed($fb, "$fo.b.bed", 1);
runCmd("intersectBed -wo -f 0.5 -a $fo.b.bed -b $fo.a.bed > $fo.o.bed");
reduce_gtb($fa, $fb, "$fo.o.bed", $fo);
runCmd("rm $fo.*.bed");

sub reduce_gtb {
  my ($fa, $fb, $fm, $fo) = @_;
  my $ta = readTable(-in => $fa, -header => 1);
  my $tb = readTable(-in => $fb, -header => 1);
  
  open(my $fhm, "<$fm") or die "cannot read $fm\n";
  my @idxsd;
  while(<$fhm>) {
    chomp;
    my ($chr1, $b1, $e1, $id, $chr2, $b2, $e2, $idx) = split "\t";
    push @idxsd, $idx;
  }
  printf "A: %5d rows, %5d removed\n", $ta->nofRow, scalar(@idxsd);
  printf "B: %5d rows added\n", $tb->nofRow;

  $ta->delRows(\@idxsd);
  $ta->rowMerge($tb);
  $ta->sort("chr", 1, 0, "beg", 0, 0);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $ta->tsv(1);
  close $fho;
}
sub gtb2bed {
  my ($fi, $fo, $opt) = @_;
  my $ti = readTable(-in => $fi, -header => 1);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  for my $i (0..$ti->lastRow) {
    my ($id, $par, $chr, $beg, $end, $srd, $locES, $locIS, $locCS, $loc5S, $loc3S, $phaseS, $src, $conf, $cat1, $cat2, $cat3, $note) = $ti->row($i);
    $cat1 eq "mRNA" || next;
    $locCS || die "no CDS for $id\n";
    
    my $rloc = [ sort {$a->[0] <=> $b->[0]} @{locStr2Ary($locCS)} ];
    my ($rtBeg, $rtEnd) = ($rloc->[0]->[0], $rloc->[-1]->[1]);
    my ($tBeg, $tEnd);
    if($srd eq "+") {
      ($tBeg, $tEnd) = ($beg+$rtBeg-1, $beg+$rtEnd-1);
    } else {
      ($tBeg, $tEnd) = ($end-$rtEnd+1, $end-$rtBeg+1);
    }
    if($opt == 1) {
      print $fho join("\t", $chr, $tBeg, $tEnd, $id)."\n";
    } elsif($opt == 2) {
      print $fho join("\t", $chr, $tBeg, $tEnd, $i)."\n";
    } else {
      die "unknown opt $opt\n";
    } 
  }
  close $fho;
}
 

__END__
