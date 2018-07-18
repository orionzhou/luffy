#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2fas.pl - convert a Gtb file to FASTA format

=head1 SYNOPSIS
  
  gtb2fas.pl [-help] [-db seq-db] [-opt option] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -d (--db)     sequence database (fasta) file
    -p (--opt)    ouput option [protein]
    -g (--usegid) use gene ID not transcript ID [No]

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use Seq;

my ($fi, $fo, $fd) = ('') x 3;
my $useGid = '';
my $opt = "protein";
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "db|d=s"   => \$fd,
  "opt|p=s"  => \$opt,
  "usegid|g"  => \$useGid,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fd;

my ($fhi, $fho);
if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot read $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my ($cnta, $cntp) = (0, 0);
my $seqH = Bio::SeqIO->new(-fh => $fho, -format => "fasta");
while( <$fhi> ) {
  chomp;
  /(^id)|(^\#)|(^\s*$)/i && next;
  my $ps = [ split("\t", $_, -1) ];
  @$ps >= 18 || die "not 19 fileds:\n$_\n";
  my ($id, $par, $chr, $beg, $end, $srd, $locES, $locIS, $locCS, $loc5S, $loc3S, $phaseS, $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
  $cat1 eq "mRNA" || next;
  $locCS || die "no CDS for $id\n";
  my $rloc = locStr2Ary($locCS);
  my $loc = $srd eq "-" ? [map {[$end-$_->[1]+1, $end-$_->[0]+1]} @$rloc] : 
    [map {[$beg+$_->[0]-1, $beg+$_->[1]-1]} @$rloc];
  $id = $par if $useGid;

  $cnta ++;
  my $seq;
  if($opt =~ /^cds$/i) {
    my $seqstr = seqRet($loc, $chr, $srd, $fd);
    $seq = Bio::Seq->new(-id=>$id, -seq=>$seqstr);
  } elsif($opt =~ /^pro/i ) {
    my $seqstr = seqRet($loc, $chr, $srd, $fd);
    my @phases = split(",", $phaseS);
    my $seqpro = Bio::Seq->new(-seq=>$seqstr)->translate(-frame=>$phases[0])->seq;
   
    $cntp ++ if $seqpro =~ /\*[\w\*]/;
    $seqpro =~ s/\*$//;
    $seqpro =~ s/\*/X/;
    die $id."\n" if !$seqpro || $seqpro =~ /^X*$/i;
    $seq = Bio::Seq->new(-id=>$id, -seq=>$seqpro);
  } elsif($opt =~ /^mrna$/) {
    $loc = [[$beg, $end]];
    my $seqstr = seqRet($loc, $chr, $srd, $fd);
    $seq = Bio::Seq->new(-id=>$id, -seq=>$seqstr);
  } elsif($opt =~ /^mrna\+$/) {
    $loc = [[$beg-1000, $end+1000]];
    my $seqstr = seqRet($loc, $chr, $srd, $fd);
    $seq = Bio::Seq->new(-id=>$id, -seq=>$seqstr);
  } else {
    die "unknown opt: $opt\n";
  }
  $seqH->write_seq($seq);
}
close $fhi;
$seqH->close();
print "$cntp / $cnta models with in-frame stop codons\n";

__END__
