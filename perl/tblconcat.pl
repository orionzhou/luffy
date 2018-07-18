#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  tblconcat.pl - concatenate multiple tabular files (with header line) into one tabular file

=head1 SYNOPSIS
  
  tblconcat.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file, can be specified multiple times
    -o (--out)    output file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;

#--------------------------------- MAIN -----------------------------------#
my @fis;
my $fo = '';
my $help_flag;
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \@fis,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !@fis || !$fo;

my $fho;
if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $n_col;
for my $i (0..$#fis) {
  my $fi = $fis[$i];
  
  open(FHI, "<$fi") or die "cannot open $fi for reading\n";
  my $header = <FHI>;
  chomp $header;
  my @colnames = split "\t", $header;

  if($i == 0) {
    print join(" ", @colnames)."\n";
    print $fho join("\t", @colnames)."\n";
    $n_col = @colnames;
  } else {
    die "inconsistent column number of $fi\n" if $n_col != @colnames;
  }

  while( <FHI> ) {
    print $fho $_;
  }
  close FHI;
}
close $fho;

exit 0;
