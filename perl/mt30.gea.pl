#!/usr/bin/perl
use strict; use Init; use Common; use Localdb; use Run; use Annotate; use Readfile;
use Bio::Seq; use Bio::SeqIO; use Path::Class; use File::Copy;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Math::Round qw/round/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
sub extractSeq {
  my ($str) = @_;
  $str =~ /\<textarea name=\"seqs\" rows=\d+ cols=\d+\>([^\<]+)/im;
  my $seqStr;
  if($1) {
    $seqStr = $1;
  } else {
    die("cannot extract seq:\n$str\n");
  }
  return $seqStr;
}
sub getProbeSeq {
  my ($fi, $fo) = rearrange(['in', 'out'], @_);
  my @ids = getId(-in=>$fi, -col=>'id');
  die scalar(@ids)." ids\n";
  my $seqH = Bio::SeqIO->new(-file=>">".$fo, -format=>"fasta");
  my $url = "http://www.plexdb.org/modules/tools/get_seq.php";
  my $cnt = 0;
  for my $id (@ids) {
    my $param = {'array_name' => 'Medicago61k', 'probesets' => $id};
    my $seqStr = extractSeq(postUrl($url, $param));
    my @tmp = split("\n", $seqStr);
    $tmp[0] =~ /^\>([\w\.\_]+): .*; ([\w\.\|\_]+)/;
    my $seq = Bio::Seq->new(-id=>$1, -seq=>$tmp[1], -description=>$2);
    $seqH->write_seq($seq);
    last if ++$cnt < 1;
  }
}
sub prepareLoad {
    my ($fIn, $fOut, $db, $table) = rearrange(["in", "out", "db", "table"], @_);
    my $fInH = new IO::File $fIn, "r";
    my $fOutH = new IO::File $fOut, "w";
    my $firstLine = readline($fInH);
    my $cntCol = split("\t", $firstLine);
    my $rst = {};
    while( <$fInH> ) {
        chomp;
        my @eleAry = split("\t", $_);
        my $probeId = $eleAry[0];
        die((@eleAry-1).".[not $cntCol] values on $_\n") unless $cntCol == @eleAry;
        my $exp = join(" ", @eleAry[1..$#eleAry]);
        print $fOutH join("\t", $probeId, $exp)."\n";
    }
    my $cmd = "LOAD DATA LOCAL INFILE '$fOut' INTO TABLE $db.$table FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';";
    print "\n$cmd\n";
}
my $DIR_Work = dir($DIR_Misc2, "mtgea");
my $f01 = file($DIR_Work, "01_seq.fa");
my $f11 = file($DIR_Work, "11_exp.txt");
my $f12 = file($DIR_Work, "12_probeset_id.txt");
#getProbeSeq(-in=>$f12, -out=>$f01);
my $f14 = file($DIR_Work, '14_db.txt');
#prepareLoad(-in=>$f11, -out=>$f14, -db=>'mt_hapmap', -table=>'expression');



