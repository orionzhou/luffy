#!/usr/bin/perl
use strict;
use Init;
use Common;
use Bio::Seq;
use Path::Class;
use Data::Dumper;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
my $DIR_Work = dir($DIR_Misc1, "chloroplast");
my $seqId = "chloroplast";
my $types = "CDS";
my $accNumGroup = [ [1..16, 19..21, 23..28, 101], [17..18, 22, 29, 30] ];
my ($feDb, $refDb) = ("mt_30") x 2;
my $ld = Localdb->new(-db=>$feDb);
my $accGroup = makeAccAry($accNumGroup);
my $accAry = [mergeArray($accGroup)];
my $rec = SeqRecover->new(-acc=>$accGroup, -covopt=>2, -fedb=>$feDb, -refdb=>$refDb);
my $loc = "chloroplast:1..124033";
my $f21 = file($DIR_Work, "21_vntS.txt");
my $f22 = file($DIR_Work, "22_vntD.txt");
#$rec->writeVntsByLoc(-loc=>$loc, -out=>$f21, -opt=>1);
$rec->writeVntsByLoc(-loc=>$loc, -out=>$f22, -opt=>2);
sub writeChloro {
    my ($locAry, $fOut, $note, $strand, $groups) = 
        rearrange(['loc', 'out', 'note', 'strand', 'groups'], @_);
    $strand ||= 1;
    my $seqid = "MtChloro";
    my $filePos = getLocIndex(-chr=>$seqid, -pos=>$locAry->[0]->[0]);
    my $locHash = { map {$seqid.":".$_->[0]."..".$_->[1] => $filePos} @$locAry };
}
my ($lenGene, $lenCds, $lenInter, $lenIntron, $lentRNA, $lenmRNA, $lenrRNA) = (0) x 7;
my @types = qw/gene mRNA tRNA rRNA CDS/;
my $locHash = {};
for my $type (@types) {
    my @feAry = getFeatures(-loc=>$seqId, -types=>$type);
    my $posAry = [ map {[$_->start, $_->end]} @feAry ];
    $locHash->{$type} = $posAry;
}
my ($oPosAry, $aPosAry, $bPosAry) = posCmp($locHash->{gene}, [[1, 124033]]);
$locHash->{intergenic} = $bPosAry;
my ($oPosAry, $aPosAry, $bPosAry) = posCmp($locHash->{mRNA}, $locHash->{CDS});
$locHash->{intron} = $aPosAry;
my $param = {
    1 => ['genome'  , [[1, 124033]]],
    2 => ['tRNA'    , $locHash->{tRNA}],
    3 => ['rRNA'    , $locHash->{rRNA}],
    4 => ['mRNA'    , $locHash->{mRNA}],
    5 => ['intergenic', $locHash->{intergenic}],
    6 => ['intron'  , $locHash->{intron}]
};
for my $key (sort {$a<=>$b} keys %$param) {
    my ($suffix, $locAry) = @{$param->{$key}};
    my $fOut = file($DIR_Work, sprintf("%02d_%s.fa", $key, $suffix));
    print join("\t", $suffix)."\n";
    #writeChloro(-loc=>$locAry, -out=>$fOut, -note=>$suffix, -groups=>$groups);
}
my @types = qw/mRNA CDS intron/;
my @mrnaAry = $ld->getFeatures(-loc=>'chloroplast', -types=>'mRNA');
for my $mrna (@mrnaAry) {
    my $mrnaId = $mrna->id;
    my $mrnaNote = join(" |", $mrna->get_tag_values("Note"));
    print "\t".join("\t", $mrna->id)."\n";
    for my $type (@types) {
        my ($posAry, $strand) = getLocStr(-id=>$mrnaId, -types=>$type);
        next if @$posAry == 0;
        my $DIR_Work2 = dir($DIR_Work, lc($type));
        system("mkdir", $DIR_Work2) unless -d $DIR_Work2;
        my $fOut = file($DIR_Work2, "$mrnaId.fa");
#writeChloro(-loc=>$posAry, -groups=>$groups, -note=>$mrnaNote, -out=>$fOut);
    }
}
