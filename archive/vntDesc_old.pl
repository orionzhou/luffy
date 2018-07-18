#!/usr/bin/perl -w
use strict; use Init; use VntEffect;
use Path::Class; use Carp::Assert; 
use constant CUTOFFCOV  => 2;
use constant CUTOFFFREQ => 0.7;
use constant CUTOFFUNIQ => 2;

my $suffix = "DEFL";
my $DIR_Work = dir($DIR_Misc, $suffix);
#mkdir($DIR_Work) unless (-d $DIR_Work);
my ($fIn, $fLoc, $fOut) = (file($DIR_Work, "vnt_snp.txt"), 
    file($DIR_In, "location_".$suffix), file($DIR_Work, "vnt_snp_effect.txt"));
print "input file doesn't exist" unless (-e $fIn && -e $fLoc);
my $fInH  = new IO::File $fIn, "r";
my $fOutH  = new IO::File $fOut, "w";
print $fOutH join("\t", "geneFamily", "geneId", "type", "region", "acc", 
    "position", "class", "effect", "desc", "aaLoc"). "\n";

my $locRef = getLoc(-in=>$fLoc, -option=>"R", -column=>"Gene");
my $cnt = 0;
if (defined $fInH) {
    while(<$fInH>) {
        chomp;
        last if ++$cnt < 0;
        my @eleAry = split("\t", $_);
        next if ($_ eq "" || $eleAry[0]=~/\D/);
        assert(@eleAry>=12);
        push (@eleAry, 0) if (@eleAry==12);
        my $vntRef = {"VntId"=>$eleAry[0], "Class"=>$eleAry[1], 
            "Region"=>$eleAry[2], "Position"=>$eleAry[3], "PosOri"=>$eleAry[4],
            "RefAllele"=>$eleAry[5], "VarAllele"=>$eleAry[6], "Accession"=>$eleAry[7],
            "NumReadsWithAllele"=>$eleAry[8], "Coverage"=>$eleAry[9],
            "Freq"=>$eleAry[10], "UniqAlns"=>$eleAry[11], "AvgQual"=>$eleAry[12]};
        if (  $vntRef->{"Coverage"}  >=  CUTOFFCOV
            &&  $vntRef->{"Freq"}      >=  CUTOFFFREQ
            &&  $vntRef->{"UniqAlns"}  >=  CUTOFFUNIQ) {
            my $infoAryRef = $locRef->{$vntRef->{"Region"}};
            my ($geneFamily, $geneId, $type) = @$infoAryRef;
            my @tmpAry = split(":", $vntRef->{"Region"});
            my $vntStr = join(":", $tmpAry[0], $vntRef->{"Position"}, 
                $vntRef->{"RefAllele"}, $vntRef->{"VarAllele"});
            my ($vntEffectRef, $vntDescRef, $aaLocRef) = vntEffect($vntStr);
            for my $i (0..scalar(@$vntDescRef)-1) {
                my ($vntEffect, $vntDesc, $aaLoc) = ($vntEffectRef->[$i], 
                    $vntDescRef->[$i], $aaLocRef->[$i]);
                print $fOutH join("\t", $geneFamily, $geneId, $type,
                    $vntRef->{"Region"}, $vntRef->{"Accession"},
                    $vntRef->{"Position"}, $vntRef->{"Class"},
                    $vntEffect, $vntDesc, $aaLoc)."\n";
            }
        }
    }
}
