#!/usr/bin/perl
use strict; use Init; use Common; use Localdb; use Run; use Annotate;
use Writefile; use Readfile; use Parser; use Vnt; use DB_File;
use Bio::TreeIO; use Bio::Seq; use Bio::SeqIO; use Path::Class;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Math::Round qw/round/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
my $co = {freq=>0.7, reads2=>2};
my $refDb = 'mt_30';
my @seqids = (0..8);
#@seqids = ('chloroplast', 'mitochondrion');
#@seqids = (0, 'tc');
for my $seqid (@seqids) {
    $seqid = "chr$seqid" if $seqid =~ /^\d$/;
    print $seqid." ...\n";
#  reformatIndel1($seqid, $refDb);
#  reformatIndel2($seqid);
#  storeVnt($seqid);
#  vntDbHead($seqid);
#  callVnt(-seqid=>$seqid, -co=>$co, -refdb=>$refDb);
#  storeCalledVnt($seqid);
#  simpleSnp2Table($seqid);
}




