#!/usr/bin/perl -w
#
#------------------------------------------------------------------------------
#                         University of Minnesota
#           Copyright 2007, Regents of the University of Minnesota
#------------------------------------------------------------------------------
# Author:
#
#  Kevin Silverstein
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN

=head1 NAME

 find_seq_motif.pl - identify positions of simple sequence motifs

=head1 SYNOPSIS

 extract_fasta.pl [-help] [-debug] [-motif regex_or_file] [-multim] 
                  [-nmm n_mismatch] [-range start:stop] [-revcomp] 
                  <fasta-file>+

  Options:
    -help       brief help message
    -motif      a single regex motif or file with multiple 1-per-line regexs
    -multim     flag to identify multiple matches in each sequence
    -nmm        number of allowed mismatches (0, 1 or 2)
    -range      restrict sequence search to a subsequence range
    -revcomp    flag to also scan for matches with the reverse comp of motifs
    -verbose    output basic status messages

=head1 DESCRIPTION

 This program identifies matches to simple regular expression patterns 
 within nucleotide-only fasta sequence files.  For now, regular expression 
 patterns are restricted to simple ones like 'ACCAG[GT][ACG]TAG'.  

=head1 OPTIONS

=over 6

=item B<-help>

 Print a usage summary.

=item B<-debug>

 Print debugging information

=item B<-motif regex_or_file>

 Either a single regular expression pattern (e.g., 'GCCA[AGT]TT[TA]CCG')
 may be specified, or alternatively a file name with one pattern per line.

 NOTE: only simple patterns consisting of residues (e.g, 'ACCGC'),
 residue options (e.g., '[AC]', '[CGT]') or combinations of the two.

 This option must be specified as there is no default pattern.

=item B<-multim>

 Enable multiple matches to be found for a single motif against each
 sequence.
 By default, only the first match is reported.

=item B<-nmm n_mismatch>

 Enable motif matches with 0, 1 or 2 residue mismatches.
 By default, no mismatches are tolerated.

=item B<-range start:stop>

 Restricts the pattern search to a subrange within each target sequence.
 The parameters start and stop are residue positions within the target
 sequence.  For example, if the following option is provided:

     -motif 'AA[GC]TCCG' -range '10:30'

 and the following sequence was provided:

 1        10        20        30        40        50
 AACTCCGTGCCACTGAAGTCCGTGTCACTACTGTCATGTGGTGCACATGCAT
 ^^^^^^^        ^^^^^^^

 Then the first motif match (pos 1-7) will be skipped, and the second
 (pos 16-22) will be reported (since only the second falls in the
 specified position range (10-30).

=item B<-revcomp>

 Flag to consider both matches to the motifs provided *and* to 
 the reverse complemented motifs as well.

=item B<fasta-file>

 Any fasta-formated file(s) or STDIN stream to use as the set of
 target sequences that are scanned for motifs.  Extracted motif
 matches and positions will be sent to STDOUT.

=back

=head1 BUGS

=head1 REFERENCES

=head1 VERSION

 0.1

=cut

#### END of POD documentation.
#-----------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;

my $help_flag;
my $debug;
my $verbose;
my $range;
my $low_pos;
my $high_pos;
my $regexp_input;
my $n_mismatch = 0;
my $multimatch_flag;
my $revcomp_flag;
my %options = (
	       "help|h"     => \$help_flag,
	       "debug|d"    => \$debug,
	       "motif=s"    => \$regexp_input,
	       "multim"     => \$multimatch_flag,
	       "nmm=s"      => \$n_mismatch,
	       "range=s"    => \$range,
	       "revcomp"    => \$revcomp_flag,
	       "verbose|v"  => \$verbose
	       );

#
# VARIABLES
#

my @regexp = ();  # Array storing each regular expression to search for
my %expanded_regexp = (); # hash for expanded n-mismatch regular expressions

#----------------------------------- MAIN -----------------------------------#
#
# Check the command-line options
#
GetOptions(%options) or pod2usage(2);
pod2usage(1) if $help_flag;

pod2usage({-message => "A simple regular expression pattern or file with multiple patterns must be specified with the -motif option.\n", -exitval => 2}) if (!$regexp_input);

if ($range) {
  pod2usage({-message => "-range option must be followed by positive integer values separated by a colon (low-pos:high-pos).\n", -exitval => 2}) if (!($range =~ /^\d+:\d+/));
  ($low_pos, $high_pos) = $range =~ /^(\d+):(\d+)/;
}

pod2usage({-message => "-nmm option must be either 0, 1 or 2.\n", -exitval => 2}) if (!($n_mismatch =~ /^[012]$/));

#
# 1. Read in the set of input motifs
#

if (-f $regexp_input) {
  &extract_regex_from_file($regexp_input, \@regexp);
}
else {
  $regexp_input =~ s/N/\[ACGT\]/g; 
  if ($regexp_input =~ /^[ACGT\[\]]+$/) {
    push @regexp, $regexp_input;
  }
  else {
    warn "WARNING: skipping unrecognized regular expression: $_\n";
  }
}

#
# Expand out the regular expression patterns to include mismatches 
#  and reverse complements (if specified).
#
&add_related_motifs($n_mismatch, $revcomp_flag, \@regexp, \%expanded_regexp);

#
# check to make sure we found at least one acceptable regexp to search for
#
if (scalar(@regexp) > 0) {
  if ($verbose || $debug) {
    print STDERR "The following motifs will be searched for:\n";
    foreach my $motif (@regexp) {
      print STDERR "  $motif";
      if ($debug && defined($expanded_regexp{$motif})) {
	foreach my $exp_motif (@{$expanded_regexp{$motif}}) {
	  print STDERR "\n      $exp_motif\n";
	}
      }
      else {
	print STDERR "\n";
      }
    }
  }
}
else {
  pod2usage({-message => "No recognizable expression patterns were found to search with.  Please carefully read instructions for the -motif option.\n", -exitval => 2});
}

#
# Systematically scan each input fasta file for each motif variant
#

print "Seq ID\tMotif\tStart Pos\tEnd Pos\tN-mismatch\tStrand\tMatch Str\n";

my $seqin; # Bio::SeqIO Fasta Input stream

if ($#ARGV < 0) { # No fasta given, expect it from STDIN
  if ($verbose || $debug) { warn "\nReading input FASTA seqs from STDIN\n";}
  $seqin = Bio::SeqIO->new(-fh => \*STDIN, -format => 'Fasta');
  &scan_motif_set($seqin, @regexp);
}
else {
  foreach my $fasta_file (@ARGV) {

    if ($verbose || $debug) { 
      warn "\nReading input FASTA seqs from $fasta_file\n";
    }

    my $fhandle;

    if ($fasta_file =~ /\.(gz|Z)$/) { # compressed file
      # Couldn't get the following to work, so I opened my own file handle!
      # $seqin = Bio::SeqIO(-file => "gzip -cd $fasta_file |", 
      #     		  -format => 'Fasta');
      open ($fhandle, "gzip -cd $fasta_file |") || 
	die "Can't uncompress $fasta_file: $!\n";
    }
    else {
      # Couldn't get the above to work, so I open my own file handles!
      # $seqin = Bio::SeqIO->new(-file => $fasta_file, -format => 'Fasta');
      open ($fhandle, "$fasta_file") || 
	die "Can't uncompress $fasta_file: $!\n";
    }

    $seqin = Bio::SeqIO->new(-fh => \*$fhandle, -format => 'Fasta');

    &scan_motif_set($seqin, @regexp);
    close(FILE);
  }
}

exit 0;

#------------------------------- SUBROUTINES --------------------------------#

sub extract_regex_from_file {
  my ($fname, $regex_arref) = @_;

  if ($verbose || $debug) { warn "\nReading list of motifs from $fname.\n";}

  open (FILE, $fname) || 
    die "Can't open file listing search motifs $fname: $!\n";
  while (<FILE>) {

    chomp;

    if (/^\s*$/) { next; }  # skip blank lines
    
    if (/^\#/) { next; } # skip comment lines

    # get rid of leading whitespace which messes up column identification
    $_ =~ s/^\s*//;

    # substitute ambiguous bases with expanded form
    s/N/\[ACGT\]/g;

    if (/^[ACGT\[\]]+$/) {
      push @{$regex_arref}, $_;
    }
    else {
      warn "WARNING: skipping unrecognized regular expression: $_\n";
    }
  }
  close(FILE);

}

sub add_related_motifs {
 my ($n_mismatch, $revcomp, $regexp_arref, $expanded_regexp_href) = @_;

 for (my $nmm=0; $nmm <=$n_mismatch; $nmm++) {
   foreach my $base_motif (@$regexp_arref) {
     my @base_motif_residues = &separate_residues($base_motif);

     my $degen_motif = &make_degen_regex($nmm, @base_motif_residues);
     push @{$$expanded_regexp_href{$base_motif}}, "$nmm\t+\t$degen_motif";

     if ($revcomp) {
       my $rc_motif = 
	 &make_degen_regex($nmm, &separate_residues(&reverse_motif(&separate_residues($base_motif))));
       push @{$$expanded_regexp_href{$base_motif}}, "$nmm\t-\t$rc_motif";
     }
   }
 }
}

sub separate_residues {
  my $motif_str = shift @_;
  my @motif_elems = ();

  # now separate the motif into terms and check them all for validity
  $motif_str =~ s/([ACGT\]])\[/$1\t\[/g;   # ATTG[ATGC]AG  => ATTG [ATGC]AG
  $motif_str =~ s/\]([\[ACGT])/\]\t$1/g;   # ATTG [ATGC]AG => ATTG [ATGC] AG
  my @subphrases = split /\t/, $motif_str;
  # if ($debug) { print STDERR "Subphrases: @subphrases\n";}  ## DEBUG
  my $term;
  foreach $term (@subphrases) {
    if ($term =~ /^\[([ACGT]+)\]$/) {
      push @motif_elems, $term;
    }
    else {
      for (my $i=0; $i<length($term); $i++) {
        my $elem = substr($term, $i, 1);
        push @motif_elems, $elem;
      }
    }
  }
  return @motif_elems;
}

sub reverse_motif {
    my @motif = @_;
    my $reversed_motif = '';

    for (my $i=$#motif; $i>=0; $i--) {
        $reversed_motif .= $motif[$i];
    }
    $reversed_motif =~ tr/acgtACGT/TGCATGCA/;

    return ($reversed_motif);
}

sub make_degen_regex {
    my ($n_mismatch, @motif_elems) = @_;

    if ($n_mismatch == 0) { return join('',@motif_elems);}

    my $motif_str = "(?:" . join('',@motif_elems);

    if ($n_mismatch == 1) { # make all combinations of 1 mismatch
        for (my $i=0; $i<=$#motif_elems; $i++) {
            my @substituted_motif = @motif_elems;
            $substituted_motif[$i] = '[ACGT]';
            $motif_str .= '|' . join('',@substituted_motif);
        }
    }
    elsif ($n_mismatch == 2) { # make all combinations of 2 mismatches
        for (my $i=0; $i<$#motif_elems; $i++) {
            my @substituted_motif = @motif_elems;
            $substituted_motif[$i] = '[ACGT]';
            for (my $j=$i+1; $j<=$#motif_elems; $j++) {
                $substituted_motif[$j] = '[ACGT]';
                $motif_str .= '|' . join('',@substituted_motif);
                $substituted_motif[$j] = $motif_elems[$j];  # set it back!
            }
        }
    }

    $motif_str .= ")";
    return $motif_str;
}

sub scan_motif_set {
  my ($seqin, @regexp) = @_;

  while (my $seqobj = $seqin->next_seq() ) {

    my $seqid = $seqobj->display_id();
    my ($seq, $init_offset);
    if ($range) {
      $seq = uc($seqobj->subseq($low_pos, $high_pos));
      $init_offset = $low_pos - 1;
    }
    else {
      $seq = uc($seqobj->seq());
      $init_offset = 0;
    }

    foreach my $base_regexp (@regexp) {
      # NOTE: as we permute through the various expanded regular expressions
      #  we only want to output a single match at a given position in each
      #  target sequence.  So we will mark each position as we see it.
      my %seen_pos = ();
      my $match_count = 0;

      foreach my $regexp_term (@{$expanded_regexp{$base_regexp}}) {
	my ($nmm, $orient, $motif) = split /\t/, $regexp_term;

	my $offset = $init_offset;
	my $rest_of_seq = $seq;

	while (($multimatch_flag && ($rest_of_seq ne '')) || 
	       (!$multimatch_flag && ($rest_of_seq ne '') &&
		($match_count == 0))) {

	  if ($debug) {warn "Scanning motif $motif ($base_regexp, $nmm, $orient) against $rest_of_seq\n";} ## DEBUG
	  if ($rest_of_seq =~ /($motif)/) {
            my $match_seq = $1;
            my $left_seg_len = length($`);
            my $low = $offset + $left_seg_len + 1;
            my $hi = $low + length($match_seq);
	    $match_count++;

	    # make sure we haven't already seen the current position
	    # earlier (i.e., with fewer mismatches)
	    if (!defined($seen_pos{$low})) {
	      print "$seqid\t$base_regexp\t$low\t$hi\t$nmm\t$orient\t$match_seq\n";
	    }

            $rest_of_seq = $match_seq . $';
	    $rest_of_seq =~ s/^\w//; # chop off first match pos

            $offset = $low; # subtract 1 to keep perl index start of 0
	    $seen_pos{$low}++;
	  }
	  else {
            $rest_of_seq = '';
	  }
	}
      }
    }
  }
}
