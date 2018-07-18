#!/usr/bin/perl -w
use strict; use Init; use Common;
use Path::Class; use Carp::Assert; 
use List::Util qw[min max sum]; use POSIX qw(ceil floor);

my $fBlat3  = file($DIR_Work, "blat_3.txt");
my ($fBlat51, $fBlat52, $fBlat53) =
    (file($DIR_Work, "blat_51.txt"), file($DIR_Work, "blat_52.txt"), file($DIR_Work, "blat_53.txt"));
#&extractPM(-in=>$fBlat3, -out1=>$fBlat51, -out2=>$fBlat52, -out3=>$fBlat53);
my ($fBlat61, $fBlat62) = (file($DIR_Work, "blat_61.txt"), file($DIR_Work, "blat_62.txt"));
my $fRawModel = file($DIR_Work, 'defl_model.txt');
#&addRawModel(-in1=>$fBlat51, -in2=>$fRawModel, -out=>$fBlat61);
my ($fBlat71, $fBlat72, $fBlat73) =
    (file($DIR_Work, "blat_71.txt"), file($DIR_Work, "blat_72.txt"), file($DIR_Work, "blat_73.txt"));
#&findGene(-in=>$fBlat61, -out1=>$fBlat71);

sub extractPM {
    my ($fIn, $fOut1, $fOut2, $fOut3) = rearrange(['in', 'out1', 'out2', 'out3'], @_);
    print "input file doesn't exist" unless (-e $fIn);
    my $fInH  = new IO::File $fIn, "r";
    my $fOutH1 = new IO::File $fOut1, "w";
    my $fOutH2 = new IO::File $fOut2, "w";
    my $fOutH3 = new IO::File $fOut3, "w";
    my ($tcHash, $pmHash, $npmHash) = ({}, {}, {});
    my ($qryId, $qryDesc);
    while(<$fInH>) {
        chomp;
        next if $_ eq "";
        my @eleAry = split("\t", $_);
        if ($eleAry[0] ne "") {
            assert @eleAry == 2;
            ($qryId, $qryDesc) = split(" ", $eleAry[1], 2);
        } else {
            assert (@eleAry == 7);
            $tcHash->{$qryId} = 0 if !exists($tcHash->{$qryId});
            $tcHash->{$qryId} ++;
            my $percentId = $eleAry[1];
            if($percentId eq "100%") {
                $pmHash->{$qryId} = 0 if ! exists $pmHash->{$qryId};
                $pmHash->{$qryId} ++;
                print $fOutH1 join("\t", $qryId, substr($qryDesc, 0))."\n" if $pmHash->{$qryId} == 1;
                print $fOutH1 "\t".join("\t", @eleAry[1..@eleAry-1])."\n";
            } else {
                $npmHash->{$qryId} = 0 if !exists $npmHash->{$qryId};
                $npmHash->{$qryId} ++;
                print $fOutH3 join("\t", $qryId, substr($qryDesc, 0))."\n" if $npmHash->{$qryId} == 1;
                print $fOutH3 "\t".join("\t", @eleAry[1..@eleAry-1])."\n";
            }
        }
    }
    my $npmHash2 = {};
    for my $acc (keys %$npmHash) {
        $npmHash2->{$acc} = 0 if !exists $pmHash->{$acc};
    }
    seek($fInH, 0, 0);
    while( <$fInH> ) {
        chomp;
        next if $_ eq "";
        my @eleAry = split("\t", $_);
        if ($eleAry[0] ne "") {
            assert @eleAry == 2;
            ($qryId, $qryDesc) = split(" ", $eleAry[1], 2);
        } else {
            assert (@eleAry == 7);
            my $percentId = $eleAry[1];
            if($percentId ne "100%" && exists($npmHash2->{$qryId})) {
                $npmHash2->{$qryId} ++;
                print $fOutH2 join("\t", $qryId, substr($qryDesc, 0))."\n" if $npmHash2->{$qryId} == 1;
                print $fOutH2 "\t".join("\t", @eleAry[1..@eleAry-1])."\n";
            }
        }
    }
    print scalar(keys %$tcHash)." TCs => ".sum(values %$tcHash)." locations, of which:\n";
    print "\t".scalar(keys %$pmHash)." TCs => ".sum(values %$pmHash)." perfect matches\n";
    print "\t".scalar(keys %$npmHash)." TCs => ".sum(values %$npmHash)." non-perfect matches\n";
    print "\t".scalar(keys %$npmHash2)." TCs => ".sum(values %$npmHash2)." non-perfect matches\n";
}
sub addRawModel {
    my ($fIn1, $fIn2, $fOut) = rearrange(["in1", "in2", "out"], @_);
    my $fInH1  = new IO::File $fIn1, "r";
    my $fOutH = new IO::File $fOut, "w";
    my $rm = &getRawModel($fIn2);
    my ($acc, $desc, $rawModel);
    my $rmHash = {};
    while( <$fInH1> ) {
        chomp;
        next if $_ eq "";
        my @eleAry = split("\t", $_);
        if ($eleAry[0] ne "") {
            assert @eleAry == 2;
            ($acc, $desc) = @eleAry[0..1];
            $rawModel = exists($rm->{$acc}) ? $rm->{$acc} : "NA";
            $rmHash->{$acc} = 0;
            if($rawModel ne "NA") {
                $rmHash->{$acc} ++;
                print $fOutH join("\t", $acc, $rawModel, $desc)."\n";
            }
        } else {
            print $fOutH $_."\n" if $rawModel ne "NA";
        }
    }
    print sum(values %$rmHash)." / ".scalar(keys %$rmHash)." TCs have supporting models\n";
}
sub findGene {
    my ($fIn, $fOut1) = rearrange(['in', 'out1'], @_);
    print "input file doesn't exist" unless (-e $fIn);
    my $fInH  = new IO::File $fIn, "r";
    my $fOutH1 = new IO::File $fOut1, "w";
    my ($cntTc, $cntLocM, $cntLoc, $cntGene, $cntGeneM) = (0, 0, 0, 0, 0);
    my ($tc2loc, $loc2gene, $geneHash, $tc2gene) = ({}, {}, {}, {});
    my ($acc, $rmPos, $desc);
    while(<$fInH>) {
        chomp;
        next if $_ eq "";
        my @eleAry = split("\t", $_);
        if ($eleAry[0] ne "") {
            assert @eleAry == 3;
            ($acc, $rmPos, $desc) = @eleAry[0..2];
            $cntTc ++;
            print $fOutH1 join("\t", @eleAry)."\n";
        } else {
            assert (@eleAry == 7);
            $cntLoc ++;
            my ($tRange) = ($eleAry[5]);
            my ($chrName, $chrStart, $chrStop) = &splitLocStr($tRange);
            print $fOutH1 "\t".join("\t", @eleAry[1..@eleAry-1])."\n";
            my @features = getFeatures($tRange, "mRNA");
            if (@features < 1) {
                print $fOutH1 "\t\t"."NA"."\n";
            } else {
                #print ">1 features: $acc\n" if (@features > 1) ;
                $cntLocM ++;
            }
            foreach my $fe (@features) {
                $cntGene ++;
                my $locStr = $fe->seq_id.":".$fe->start."..".$fe->end;
                $tc2loc->{$acc} = [] if !exists($tc2loc->{$acc});
                push (@{$tc2loc->{$acc}}, $locStr);
                $loc2gene->{$locStr} = [] if !exists($loc2gene->{$locStr});
                push (@{$loc2gene->{$locStr}}, $fe->display_name);
                my $fePosStr = &getFePos($fe);
                my ($loc1, $loc2) = ($chrStart, $chrStop);
                $loc1 = $fe->start if $loc1 < $fe->start;
                $loc2 = $fe->stop if $loc2 > $fe->stop;
                my $overlap = $loc2 - $loc1 + 1;
                my $lenStr = ($chrStop-$chrStart+1)."/".$overlap."/".$fe->length;
                $geneHash->{$fe->id} = 0 if ! exists($geneHash->{$fe->id});
                $geneHash->{$fe->id} ++;
                if(&cmpPos(-parent=>$fePosStr, -child=>$rmPos) == 1) {
                    $tc2gene->{$acc} = [] if !exists($tc2gene->{$acc});
                    push (@{$tc2gene->{$acc}}, $fe->id);
                    $cntGeneM ++;
                }
                my $note = getNote($fe->id);
                print $fOutH1 "\t\t".join("\t", $fePosStr, $lenStr, $fe->display_name, $locStr, $note)."\n";
            }
        }
    }
    assert($cntGene == sum(values %$geneHash));
    print $cntTc." TCs -> ".$cntLoc." locations, of which: \n";
    print "\t".scalar(keys %$tc2loc)." TCs / ".$cntLocM."(".scalar(keys %$loc2gene)." unique)"
        ." locations -> ".$cntGene."(".scalar(keys %$geneHash)." unique) IMGAG genes\n";
    #print "\t".scalar(keys %$tc2gene)." TCs are correctly mapped to ".$cntGeneM." genes\n";
}



sub getRawModel {
    my ($fIn) = @_;
    my $fInH = new IO::File $fIn, "r" || die("cannot open file: $fIn\n");
    my $rst = {};
    while( <$fInH> ) {
        chomp;
        next if !$_;
        my @eleAry = split("\t", $_);
        my $acc = $eleAry[0];
        my @posAry;
        for my $ele (@eleAry[1..$#eleAry]) {
            push (@posAry, $ele) if $ele;
        }
        my $flag = 0;
        $flag = 1 if @posAry % 2 != 0;
        my @posStrAry;
        for my $i (0..$#posAry) {
            $flag = 1 if $i %2 == 0 && $i > 1 && $posAry[$i] - $posAry[$i-1] != 1;
            push (@posStrAry, join("..", $posAry[$i-1], $posAry[$i])) if $i % 4 == 1; #caution!
        }
        $rst->{$acc.".1"} = join(";", @posStrAry);
        if($flag) {
            print "Error!\n".$_."\n";
            exit 1;
        }
    }
    return $rst;
}
