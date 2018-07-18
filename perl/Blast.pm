package Blast;
use strict; 
use Common;
use Location;
use Seq;
use Bio::SeqIO;
use File::Basename;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use URI::Escape;
use LWP::UserAgent;
use HTTP::Request::Common qw/POST/;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/blast2psl blast2gal
  parse_aln_string
  split_blast_output blast_tiling blast_nr blast_nr_batch/;
@EXPORT_OK = qw//;

sub parse_aln_string {
  my ($qSeq, $tSeq) = @_;
  my (@qPos, @tPos, @qIns, @tIns);
  my ($qLen, $tLen) = (length($qSeq), length($tSeq));
  $qLen == $tLen || die "unequal seqlen: $qLen <> $tLen\n$qSeq\n";

  my ($qPos, $tPos) = (1, 1);
  for my $i (1..$qLen) {
    my ($chQ, $chT) = (substr($qSeq, $i-1, 1), substr($tSeq, $i-1, 1));
    push @qPos, $qPos;
    push @tPos, $tPos;
    if($chQ =~ /[\-\_\s]/) {
      push @tIns, [$i, $i];
    } else {
      $qPos ++;
    }
    if($chT =~ /[\-\_\s]/) {
      push @qIns, [$i, $i];
    } else {
      $tPos ++;
    }
  }
  my ($qIns, $tIns) = (scalar(@qIns), scalar(@tIns));
  my ($qInsL, $tInsL) = (posMerge(\@qIns), posMerge(\@tIns));
  my $qNumIns = $qInsL ? scalar(@$qInsL) : 0;
  my $tNumIns = $tInsL ? scalar(@$tInsL) : 0;
  my ($uLoc) = posDiff([[1, $qLen]], [@qIns, @tIns]);

  my (@qLoc, @tLoc, @stats);
  for (@$uLoc) {
    my ($idxB, $idxE) = @$_;
    push @qLoc, [$qPos[$idxB-1], $qPos[$idxE-1]];
    push @tLoc, [$tPos[$idxB-1], $tPos[$idxE-1]];
    my $qSeq1 = substr($qSeq, $idxB-1, $idxE-$idxB+1);
    my $tSeq1 = substr($tSeq, $idxB-1, $idxE-$idxB+1);
    die "$qSeq1\n" if $qSeq1 =~ /[\-\_]/;
    die "$tSeq1\n" if $tSeq1 =~ /[\-\_]/;
    push @stats, [seqCompare($qSeq1, $tSeq1)];
  }
  return (\@qLoc, \@tLoc, \@stats, $qNumIns, $qIns, $tNumIns, $tIns);
}
sub blast2psl {
  my ($ps) = @_;
  my ($qId, $qBeg, $qEnd, $qSize, $tId, $tBeg, $tEnd, $tSize, $alnLen, 
    $match, $misMatch, $gaps, $e, $score, $qSeq, $tSeq) = @$ps;
  $qBeg <= $qEnd || die "$qId $qBeg > $qEnd\n";
  $alnLen == $match+$misMatch+$gaps || die "len error\n".join("\t", @$ps)."\n";
  my ($qSrd, $tSrd) = ("+") x 2;
  if($tBeg > $tEnd) {
    ($tBeg, $tEnd) = ($tEnd, $tBeg);
    $qSrd = "-";
  }

  my ($qLoc, $tLoc, $stat, $qNumIns, $qIns, $tNumIns, $tIns) = parse_aln_string($qSeq, $tSeq);
#  print join("\t", $qId, $qBeg, $qEnd, $qSize, $tId, $tBeg, $tEnd, $tSize)."\n";
  my $nBlock = @$qLoc;
  my (@qBegs, @tBegs, @blockLens);
  for my $i (0..$nBlock-1) {
    my ($qb, $qe) = @{$qLoc->[$i]};
    my ($tb, $te) = @{$tLoc->[$i]};
    my $qLen = $qe - $qb + 1;
    my $tLen = $te - $tb + 1;
    $qLen == $tLen || die "len error: $qb-$qe $tb-$te\n";

    my ($tBegF, $qBegF);
    if($qSrd eq "-") {
      $tBegF = $tEnd - $te + 1;
      $qBegF = $qSize-($qBeg+$qe-1)+1;
    } else {
      $tBegF = $tBeg + $tb - 1;
      $qBegF = $qBeg + $qb - 1;
    }
    push @blockLens, $qLen;
    push @tBegs, $tBegF - 1;
    push @qBegs, $qBegF - 1;
#    print join("\t", $qb, $qe, $tb, $te, $qLen)."\n";
#    print join("\t", $qBegF, $tBegF)."\n";
#    die if $i == 2;
  }
  my $qBegStr = join(",", @qBegs);
  my $tBegStr = join(",", @tBegs);
  my $blockLenStr = join(",", @blockLens);
  
  return [ $match, $misMatch, 0, 0,
    $qNumIns, $qIns, $tNumIns, $tIns, $qSrd,
    $qId, $qSize, $qBeg-1, $qEnd, $tId, $tSize, $tBeg-1, $tEnd, 
    $nBlock, $blockLenStr, $qBegStr, $tBegStr ];
}


sub write_blast {
  my ($id, $rows, $fo) = @_;
  open(my $fho, ">$fo") or die "cannot open $fo for writing\n";
  print $fho join("\t", qw/qId qBeg qEnd strand hId hBeg hEnd pct e score/)."\n";
  print $fho join("\n", map {join("\t", @$_)} @$rows)."\n";
  close $fho;
}
sub split_blast_output {
  my ($fi, $dirO) = @_;
  make_path($dirO) unless -d $dirO;
  
  open(my $fhi, "<$fi") or die "cannot open $fi for reading\n";
  my ($idP, $lineP) = ("", "");
  my @rows;
  while(my $line = <$fhi>) {
    chop($line);
    my ($qId, $hId, $pct, $alnLen, $mm, $gap, $qBeg, $qEnd, $hBeg, $hEnd, $e, $score) = split("\t", $line);
    my $strand = $hBeg > $hEnd ? "-" : "+";
    ($hBeg, $hEnd) = ($hEnd, $hBeg) if $strand eq "-";
    my $row = [$qId, $qBeg, $qEnd, $strand, $hId, $hBeg, $hEnd, $pct, $e, $score];
    if($qId eq $idP) {
      push @rows, $row;
    } else {
      write_blast($idP, \@rows, "$dirO/$idP.tbl") if $idP;
      @rows = ($row);
      $idP = $qId;
    }
  }
  write_blast($idP, \@rows, "$dirO/$idP.tbl") if @rows > 0;
}


sub blast_tiling {
  my ($fl, $dirI, $fo) = @_;
  open(my $fho, ">$fo") or die "cannot open $fo for writing\n";
  print $fho join("\t", qw/qId qBeg qEnd strand hId hBeg hEnd qLen hLen pct e score/)."\n";
 
  my $tl = readTable(-in=>$fl, -header=>1);
  my @ids = $tl->col("id");
  my ($cntB, $cntG) = (0, 0);
  for my $id (sort @ids) {
    my $fi = "$dirI/$id.tbl";
    if( ! -s $fi ) {
      $cntB ++;
      next;
    }
   
    my $ti = readTable(-in=>$fi, -header=>1);
    my @locs = map { [$ti->elm($_, "qBeg"), $ti->elm($_, "qEnd")] } (0..$ti->nofRow-1);
    my @es = $ti->col("e");
    my $refs = tiling(\@locs, \@es, 1);
    for (@$refs) {
      my ($beg, $end, $idx) = @$_;
      next if ($end - $beg + 1) < 100;
      my ($qId, $qBeg, $qEnd, $strd, $hId, $hBeg, $hEnd, $pct, $e, $score) = $ti->row($idx);
      my $begH = sprintf "%d", $hBeg + ($beg-$qBeg) * ($hEnd-$hBeg)/($qEnd-$qBeg);
      my $endH = sprintf "%d", $hEnd - ($qEnd-$end) * ($hEnd-$hBeg)/($qEnd-$qBeg);
      my $qLen = $end - $beg + 1;
      my $hLen = $endH - $begH + 1;
      print $fho join("\t", $id, $beg, $end, $strd, $hId, $begH, $endH, $qLen, $hLen, $pct, $e, $score)."\n";
    }
    printf " %4d: %10d\n", (++$cntG)+$cntB, $ti->nofRow;
  }
  close $fho;
  printf "\n %4d / %4d tiled\n", $cntG, $cntG+$cntB;
}
sub blast_nr {
  my ($id, $seq, $fo) = @_;
  my $ua = LWP::UserAgent->new;

  my $qry = uri_escape(">$id\n$seq\n");

  my ($program, $db) = ("blastn", "nr");
  my $args = "CMD=Put&PROGRAM=$program&DATABASE=$db&QUERY=".$qry;
  my $req = new HTTP::Request POST => 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi';
  $req->content_type('application/x-www-form-urlencoded');
  $req->content($args);

  my $response = $ua->request($req);
  $response->content =~ /^    RID = (.*$)/m;
  my $rid=$1;
  $response->content =~ /^    RTOE = (.*$)/m;
  my $rtoe=$1;
  printf "  %s: RID [%s] [%5ds]\n", $id, $rid, $rtoe;
  sleep $rtoe;

  while(1) {
      sleep 3;

      $req = new HTTP::Request GET => "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";
      $response = $ua->request($req);

      if ($response->content =~ /\s+Status=WAITING/m) {
          # print STDERR "Searching...\n";
          next;
      }
      if ($response->content =~ /\s+Status=FAILED/m) {
          print STDERR "Search $rid failed; please report to blast-help\@ncbi.nlm.nih.gov.\n";
          exit 4;
      }
      if ($response->content =~ /\s+Status=UNKNOWN/m) {
          print STDERR "Search $rid expired.\n";
          exit 3;
      }
      if ($response->content =~ /\s+Status=READY/m) {
          if ($response->content =~ /\s+ThereAreHits=yes/m) {
              #  print STDERR "Search complete, retrieving results...\n";
              last;
          } else {
              print STDERR "No hits found.\n";
              exit 2;
          }
      }

      # if we get here, something unexpected happened.
      exit 5;
  }

  $req = new HTTP::Request GET => "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=$rid";
  $response = $ua->request($req);

  open (FHO, ">$fo") || die "Can't open file $fo for writing: $!\n";
  print FHO $response->content;
  close FHO;
}
sub blast_nr_batch {
  my ($fi, $do) = @_;
  make_path($do) unless -d $do;

  my $seqH = Bio::SeqIO->new(-file=>"<$fi", -format=>'fasta');
  while(my $seqO = $seqH->next_seq()) {
    my ($id, $seq) = ($seqO->id, $seqO->seq);
    my $fo = "$do/$id.txt";
    if( -s $fo ) {
      print "  $id: exists - skipped\n";
      next;
    }
    blast_nr($id, $seq, $fo);
  }
  $seqH->close();
}

1;
__END__
