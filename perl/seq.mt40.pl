#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seq.mt40.pl - 

=head1 SYNOPSIS
  
  seq.mt40.pl [-help] 

  Options:
    -h (--help)   brief help message
    -i (--in)     input sequence file (fasta)
    -o (--out)    output sequence file (fasta)

=cut
  
#### END of POD documentation.
#------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::SeqIO;
use Gal;

my ($fi, $fo, $fp) = ('') x 3;
my $gap = 1000;
my $help_flag;

#------------------------------ MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "path|p=s" => \$fp,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fp;

open(my $fhi, "<$fi") || die "cannot read $fi";
open(my $fho, ">$fo") || die "cannot write $fo";
open(my $fhp, ">$fp") || die "cannot write $fp";

my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');
my $hc = {'chrU' => "scaffold"};

my $h;
while(my $seqO = $seqHI->next_seq()) {
  my ($id, $len, $seq) = ($seqO->id, $seqO->length, $seqO->seq);

  my $flag = 0;
  for my $nid (keys(%$hc)) {
    my $pre = $hc->{$nid};
    if($id =~ /^$pre/) {
      $h->{$nid} ||= [];
      my $beg = @{$h->{$nid}} ? $h->{$nid}->[-1]->[1] + $gap + 1 : 1;
      my $end = $beg + $len - 1; 
      push @{$h->{$nid}}, [$beg, $end, $id, $len];
      $flag = 1;
    }
  }
  if($flag == 0) {
    $h->{$id} = 0;
  }
}
$seqHI->close();

my $db = Bio::DB::Fasta->new($fi);
my $seqHO = Bio::SeqIO->new(-fh=>$fho, -format=>'fasta');
print $fhp join("\t", @HEAD_GAL)."\n";
my $i = 0;
for my $id (sort(keys(%$h))) {
  if($h->{$id} == 0) {
    my $len = $db->length($id);
    print $fhp join("\t", $i++, $id, 1, $len, "+", $len,
      $id, 1, $len, "+", $len, 
      '', $len, $len, 0, 0, 0, '', '', "1-$len", "1-$len")."\n";
    $seqHO->write_seq(Bio::Seq->new(-id=>$id, -seq=>$db->seq($id)));
  } else {
    my $seqstr = "";
    my $size = $h->{$id}->[-1]->[1];
    for (@{$h->{$id}}) {
      my ($beg, $end, $sid, $len) = @$_;
      print $fhp join("\t", $i++, $id, $beg, $end, "+", $size,
        $sid, 1, $len, "+", $len, 
        '', $len, $len, 0, 0, 0, '', '', "$beg-$end", "1-$len")."\n";
      $seqstr .= "N" x $gap if $seqstr ne "";
      $seqstr .= $db->seq($sid);
    }
    $seqHO->write_seq(Bio::Seq->new(-id=>$id, -seq=>$seqstr));
  }
}
$seqHO->close();
close $fhp;

exit 0;


