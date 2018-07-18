#!/usr/bin/perl -w
=pod BEGIN
  
=head1 NAME
  
  gffemsembl.pl: fix GFF3 files from Ensembl

=head1 SYNOPSIS
  
 gffemsembl.pl [-help] -i <input-file> -o <output-file>

 Options:
   -h (--help)    brief help message
   -i (--in)      input file
   -o (--out)     output file

=cut

use strict;
use FindBin;
use lib "$FindBin::Bin";

use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw/make_path/;
use Common;
use Gff;

my $help_flag;
my ($fi, $fo) = ('', '');
my %options = (
  "help|h" => \$help_flag,
  "in|i=s" => \$fi,
  "out|o=s" => \$fo,
);

#----------------------------------- MAIN -----------------------------------#
GetOptions(%options) or pod2usage(2);
pod2usage(1) if $help_flag;

my ($fhi, $fho);
if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot readn $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $hi;
my $hf;
while(<$fhi>) {
  chomp;
  if(/^\#/) {
    print $fho $_."\n";
  } else {
    my @ps = split "\t";
    my ($chr, $src, $type, $beg, $end, $score, $srd, $phase, $tagstr) = @ps;
    my $ht = parse_gff_tags($tagstr);
    if($type eq "protein_coding_gene") {
      $type = "gene";
      $hi->{$ht->{"ID"}} = 'mRNA';
    } elsif($type eq "rRNA_gene") {
      $type = "gene";
      $hi->{$ht->{"ID"}} = "rRNA";
    } elsif($type eq "tRNA_gene") {
      $type = "gene";
      $hi->{$ht->{"ID"}} = "tRNA";
    }
    if(($type eq "exon" || $type eq "CDS") && exists $ht->{"ID"}) {
      delete $ht->{"ID"};
    }
    
    if(exists $ht->{'ID'} && exists $ht->{'Parent'} && $ht->{'ID'} eq $ht->{'Parent'}) {
      my $oid = $ht->{"ID"};
      my $nid = "$oid.1";
      $hf->{$oid} = $nid;
      $ht->{"ID"} = $nid;
    }
    if($type eq "exon" && exists $ht->{"Parent"} && $hf->{$ht->{"Parent"}}) {
      $ht->{"Parent"} = $hf->{$ht->{"Parent"}};
    }
    if(exists $ht->{"Parent"} && exists $hi->{$ht->{"Parent"}}) {
      $type = $hi->{$ht->{"Parent"}};
    }
    my $ntagstr = join(";", map {$_."=".$ht->{$_}} keys(%$ht));
    @ps[2,8] = ($type, $ntagstr);
    print $fho join("\t", @ps)."\n";
  }
}
close $fhi;
close $fho;


