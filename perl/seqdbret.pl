#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqdbret.pl - retrieve given sequences (with specified ranges) from an indexed fasta file 

=head1 SYNOPSIS
  
  seqdbret.pl [-help] [-db sequence-db] [-win window-file] [-out output]

  Options:
    -h (--help)   brief help message
    -d (--db)     sequence db
    -w (--win)    window file containing ranges (can be 'stdin')
    -o (--out)    output (can be 'stdout')

=head1 BUGS
  
=head1 REFERENCES
  
=head1 VERSION
  
  0.1
  
=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;

my ($fd, $fw, $fo) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "db|d=s"   => \$fd,
  "win|w=s"  => \$fw,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fd;
my ($fhw, $fho);

if ($fw eq "" || $fw eq "stdin" || $fw eq "-") {
  $fhw = \*STDIN;
} else {
  open ($fhw, $fw) || die "Can't open file $fw: $!\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $db = Bio::DB::Fasta->new($fd);    
my $seqHO = Bio::SeqIO->new(-fh=>$fho, -format=>'fasta');

while( <$fhw> ) {
  chomp;
  next if /(^#)|(^id\s)|(^chr\s)/;
  my @ps = split "\t";
  next unless @ps >= 1;
  my ($chr, $beg, $end) = ($ps[0], '', '');
  if(@ps == 1) {
    ($beg, $end) = (1, $db->length($chr));
  } elsif(@ps == 2) {
    ($beg, $end) = (1, $ps[1]);
  } else {
    ($beg, $end) = @ps[1,2];
  }
  my $seqStr = $db->seq($chr, $beg, $end);
  $seqHO->write_seq( Bio::Seq->new(-id=>"$chr:$beg-$end", -seq=>$seqStr) );
}
$seqHO->close();

exit 0;



