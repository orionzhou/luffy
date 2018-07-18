package VntOut;

=head1 NAME

=cut
use Common; use strict; use Switch;
use Path::Class; use DBI; use Carp::Assert; use IO::File;
use List::Util qw[min max sum]; use POSIX qw(ceil floor);

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw(vntOut);
@EXPORT_OK = qw();

=head2 new
=cut
sub new {
    my $self = shift;
    my ($fIn, $fLoc, $cutOffCov, $cutOffFreq, $cutOffUniq, $cutOffMAC, $cutOffGt, $windowSize, $offSet) =
        rearrange(["in", "fLoc", "cutOffCov", "cutOffFreq", "cutOffUniq", "cutOffMAC", "cutOffGt", "windowSize", "offSet"], @_);
    $cutOffCov = $cutOffCov ? $cutOffCov:2;
    $cutOffFreq = $cutOffFreq ? $cutOffFreq:0.7;
    $cutOffUniq = $cutOffUniq ? $cutOffUniq:2;
    $cutOffMAC = $cutOffMAC ? $cutOffMAC:1;
    $windowSize = $windowSize ? $windowSize:-1;
    $offSet = $offSet ? $offSet:-1;
    my $fInH = new IO::File $fIn, "r";
    my $fLocH = new IO::File $fLoc, "w";
    my $accAry = &getAcc($fInH);
    return bless { cutOffCov =>$cutOffCov,
                      cutOffFreq=>$cutOffFreq,
                      cutOffUniq=>$cutOffUniq,
                      cutOffMAC =>$cutOffMAC,
                      cutOffGt  =>$cutOffGt,
                      windowSize=>$windowSize,
                      offSet    =>$offSet,
                      accAry    =>$accAry,
                      inhandle  =>$fInH,
                      lochandle =>$fLocH,
                      cntVnt    =>0,
                      cntMono   =>0,
                      cntBi     =>0,
                      cntMulti  =>0,
                      cntPicked =>0,
                      posCntRef =>{},
                      accSeqRef =>{},
                      posPickedRef=>[],
                      round     =>0},
        ref($self) || $self;
}

=item
=cut
sub vntOut {
    my $self = shift;
    my ($outFormat, $fOut) = @_;
    my $accAry = $self->{accAry};
    my $posPrev = 0;
    my %accNt;
    my ($startNext) = (0);
    while(1) {
        my $vntAryRef = $self->getOnePos();
        if($vntAryRef) {
            assert(scalar(@$vntAryRef) % scalar(@$accAry) == 0);
            my @tmp = split("\t", $vntAryRef->[0]);
            my $posNow = $tmp[3];
            assert($posPrev != $posNow);
        
            %accNt = $self->vntExtract($vntAryRef, $posNow);
            my ($alleleCnt, $MAC, $percentGt) = &vntFilter(%accNt);
            $self->statUpdate($alleleCnt, $MAC, $percentGt, \%accNt, $posNow);
            
            my ($cntVnt, $cntMono, $cntBi, $cntMulti, $cntPicked) = 
                ($self->{cntVnt}, $self->{cntMono}, $self->{cntBi}, $self->{cntMulti}, $self->{cntPicked});
            my $cntPos = scalar(keys %{$self->{posCntRef}});
            #print join("*", sum($cntMono, $cntBi, $cntMulti), $cntPos, $cntVnt, sum(values %{$self->{posCntRef}})), "\n";
            
            if($cntPicked == $self->{offSet} && $startNext == 0) {
                $startNext = tell($self->{inhandle});
                #print "--".$startNext."\n";
            }
            if($cntPicked == $self->{windowSize}) {
                my $offSet = $self->{offSet}>0 && $self->{offSet}<$cntPicked ? $self->{offSet}:$cntPicked;
                #print join("|", sum($cntMono, $cntBi, $cntMulti), $cntPos, $cntVnt, $cntPicked)."\n";
                $self->vntSum;
                $self->vntWrite($outFormat, $fOut);
                $self->vntClear;
                assert($startNext != 0);
                seek($self->{inhandle}, $startNext-tell($self->{inhandle}), 1);
                #print "--".tell($self->{inhandle})."\n";
                $startNext = 0;
            }
            $posPrev = $posNow;
        } else {
            $self->vntSum;
            $self->vntWrite($outFormat, $fOut);
            $self->vntClear;
            last;
        }
        last if($self->{cntVnt} < 0);
    }
}
sub vntSum {
    my $self = shift;
    my $fLocH = $self->{lochandle};
    my ($cntVnt, $cntMono, $cntBi, $cntMulti, $cntPicked, $posCntRef, $accSeqRef, $posPickedRef) = 
        ($self->{cntVnt}, $self->{cntMono}, $self->{cntBi}, $self->{cntMulti}, $self->{cntPicked},
          $self->{posCntRef}, $self->{accSeqRef}, $self->{posPickedRef});
    $self->{round} ++;
    my $round = $self->{round};
    my $cntPos = scalar(keys %$posCntRef);
    #print sum($cntMono, $cntBi, $cntMulti)."/".$cntPos." ".$cntVnt."/".sum(values %$posCntRef)."\n";
    #print $cntPicked."/".scalar(@$posPickedRef)."\n";
    assert($cntVnt == sum(values %$posCntRef) && $cntPicked == scalar(@$posPickedRef));
    my ($snpStart, $snpStop) = (($self->{offSet})*($round-1)+1, ($self->{offSet})*($round-1)+$self->{windowSize});
    $snpStop = ($self->{windowSize}==$cntPicked) ? $snpStop : ($snpStart+$cntPicked-1);
    my ($posStart, $posStop) = (min(keys %$posCntRef), max(keys %$posCntRef));
    print $fLocH join("\t", $round, $snpStart, $snpStop, $posStart, $posStop, $snpStop-$snpStart+1, $posStop-$posStart+1), "\n";
    print "round ".sprintf("%03d", $round)." : $snpStart - $snpStop [$posStart - $posStop]\n";
    print $cntVnt." Variants / ".$cntPos." UniqLoci [";
    print $cntMono."(mono-allelic) + ".$cntBi."(bi-allelic) + ".$cntMulti."(multi-allelic)]\n";
    print $cntPicked." picked with cutOff(Minor Allele Count)=".$self->{cutOffMAC}."\n";
    print "================================\n";
}
sub vntClear {
    my $self = shift;
    ($self->{cntVnt}, $self->{cntMono}, $self->{cntBi}, $self->{cntMulti}, $self->{cntPicked},
          $self->{posCntRef}, $self->{accSeqRef}, $self->{posPickedRef}) =
        (0, 0, 0, 0, 0, {}, {}, []);
}
sub statUpdate {
    my $self = shift;
    my ($alleleCnt, $MAC, $percentGt, $accNtRef, $posNow) = @_;
    if($alleleCnt <= 1) {
        $self->{cntMono} ++;
    } elsif($alleleCnt == 2) {
        $self->{cntBi} ++;
        if($MAC>=$self->{cutOffMAC} && $percentGt >=$self->{cutOffGt}) {
            for my $accTmp (keys %$accNtRef) {
                $self->{accSeqRef}->{$accTmp} .= $accNtRef->{$accTmp};
            }
            push (@{$self->{posPickedRef}}, $posNow);  #
            $self->{cntPicked} ++;
        }
    } else {
        $self->{cntMulti} ++;
    }
}
sub vntExtract {
    my $self = shift;
    my ($vntAryRef, $posNow) = @_;
    my %accNt;
    my $vntInPos = scalar(@$vntAryRef) / scalar(@{$self->{accAry}});
    $self->{posCntRef}->{$posNow} = $vntInPos; #update %posCntRef $cntVnt
    $self->{cntVnt} += $vntInPos;
    for my $vntRow (@$vntAryRef) {
        my @eleAry = split("\t", $vntRow);
        assert(@eleAry >= 12);
        push (@eleAry, 0) if (@eleAry==12);
        my $vntRef = {"VntId"=>$eleAry[0], "Class"=>$eleAry[1], 
                "Region"=>$eleAry[2], "Position"=>$eleAry[3], "PosOri"=>$eleAry[4],
                "RefAllele"=>$eleAry[5], "VarAllele"=>$eleAry[6], "Accession"=>$eleAry[7],
                "NumReadsWithAllele"=>$eleAry[8], "Coverage"=>$eleAry[9],
                "Freq"=>$eleAry[10], "UniqAlns"=>$eleAry[11], "AvgQual"=>$eleAry[12]};
        my $acc = $vntRef->{"Accession"};
        assert($posNow == $vntRef->{"Position"});
        if ( $vntRef->{"Coverage"} >= $self->{cutOffCov} ) {
            if ($vntRef->{"Freq"} >= $self->{cutOffFreq} &&  $vntRef->{"UniqAlns"} >= $self->{cutOffUniq}) {
                if(!exists($accNt{$acc})) {
                    $accNt{$acc} = uc($vntRef->{"VarAllele"});
                } elsif($accNt{$acc} eq $vntRef->{"RefAllele"}) {
                    $accNt{$acc} = uc($vntRef->{"VarAllele"});
                } else {
                    print "\tHetero: pos=".$vntRef->{"Position"}." acc=".$acc." ".
                        $vntRef->{"RefAllele"}."->".$vntRef->{"VarAllele"}."|".$accNt{$acc}."\n";
                }
            } else {
                $accNt{$acc} = uc($vntRef->{"RefAllele"}) if(!exists($accNt{$acc}));
                assert($accNt{$acc} ne "N");
            }
        } else {
            $accNt{$acc} = "N" if(!exists($accNt{$acc}));
            assert($accNt{$acc} eq "N");
        }
    }
    return %accNt;
}
sub vntFilter {
    my (%accNt) = @_;
    my %vntHash;
    for my $acc (keys %accNt) {
        my $allele = $accNt{$acc};
        if($allele ne "N") {
            $vntHash{$allele} = 0 if (!exists($vntHash{$allele}));
            $vntHash{$allele} ++;
        }
    }
    my $alleleCnt = scalar(keys %vntHash);
    my ($MAC, $percentGt);
    if($alleleCnt == 0) {
        ($MAC, $percentGt) = (0, 0);
    } else {
        $MAC = min(values %vntHash);
        $percentGt = sum(values %vntHash) / scalar(keys %accNt);  
    }
    return ($alleleCnt, $MAC, $percentGt);
}
sub getOnePos {
    my $self = shift;
    my $fInH = $self->{inhandle};
    my $accAry = $self->{accAry};
    my @vntLines;
    my $cnt = 0;
    my ($lineLen, $tag) = (0, 0);
    assert(defined $fInH);
    my $pos;
    while(my $line = readline($fInH)) {
        $lineLen = length($line);
        chomp($line);
        next if($line eq "");
        my @eleAry = split("\t", $line);
        assert(@eleAry>=12);
        $cnt ++;
        if($cnt == 1) {
            $pos = $eleAry[3];
            push (@vntLines, $line);
        } elsif($eleAry[3] == $pos) {
            push (@vntLines, $line);
        } else {
            $tag = 1;
            last;
        }
    }
    if($tag == 1) {
        seek($fInH, -$lineLen, 1);
    }
    if(scalar(@vntLines) / scalar(@{$self->{accAry}}) >= 1) {
        return \@vntLines;
    } else {
        return 0;
    }
}
sub getAcc {
    my ($fInH) = @_;
    my @accAry;
    my $cnt = 0;
    assert(defined $fInH);
    my $firstline = <$fInH>;
    my ($firstVnt) = ("");
    while(<$fInH>) {
        chomp;
        my @eleAry = split("\t", $_);
        next if ($_ eq "" || $eleAry[0]=~/\D/);
        assert(@eleAry>=12);
        push (@eleAry, 0) if (@eleAry==12);
        my $vntRef = {"VntId"=>$eleAry[0], "Class"=>$eleAry[1], 
            "Region"=>$eleAry[2], "Position"=>$eleAry[3], "PosOri"=>$eleAry[4],
            "RefAllele"=>$eleAry[5], "VarAllele"=>$eleAry[6], "Accession"=>$eleAry[7],
            "NumReadsWithAllele"=>$eleAry[8], "Coverage"=>$eleAry[9],
            "Freq"=>$eleAry[10], "UniqAlns"=>$eleAry[11], "AvgQual"=>$eleAry[12]};
        if($cnt == 0) {
            ($firstVnt) = ($vntRef->{"VntId"});
            push (@accAry, $vntRef->{"Accession"});
        } else {
            if($vntRef->{"VntId"} == $firstVnt) {
                push (@accAry, $vntRef->{"Accession"});
            } else {
                last;
            }
        }
        $cnt ++;
    }
    assert($cnt == @accAry);
    print $cnt." accessions:\n";
    print join("\t", @accAry)."\n";
    print "================================\n";
    seek($fInH, 0, 0);
    $firstline = <$fInH>;
    return \@accAry;
}
sub vntWrite {
    my $self = shift;
    my ($outFormat, $fPrefix) = @_;
    my $round = $self->{round};
    my $fOut = sprintf("%s_%03d", $fPrefix, $round);
    if($outFormat == 1) {
        $self->vntWrite_haploview($fOut);
    } elsif($outFormat == 2) { 
        $self->vntWrite_phase($fOut);
    } elsif($outFormat == 3) { 
        $self->vntWrite_rsq($fOut);
    } elsif($outFormat == 4) { 
        $self->vntWrite_LDhat($fOut);
    } elsif($outFormat == 5) { 
        $self->vntWrite_LDhot($fOut);
    } else {
        print "error outFormat".$outFormat."\n"; exit 1;
    }
}
sub vntWrite_haploview {
    my $self = shift;
    my ($fOut) = @_;
    my ($posAryRef, $accSeqRef) = ($self->{posPickedRef}, $self->{accSeqRef});
    my $fOutH1 = new IO::File $fOut.".txt", "w";
    my $fOutH2 = new IO::File $fOut."_loc.txt", "w";
    my $cntPos = scalar(@$posAryRef);
    my $accAryRef = $self->{accAry};
    my $cntAcc = scalar(@$accAryRef);
    my $eHash1 = {"A"=>1, "C"=>2, "G"=>3, "T"=>4, "N"=>0};
    my $eHash2 = {"A"=>"A", "C"=>"C", "G"=>"G", "T"=>"T", "N"=>"?"};
    my $eHash3 = {"A"=>"L", "C"=>"O", "G"=>"P", "T"=>"Q", "N"=>"Z"};
    my $ttt = 0;
    for my $acc (@$accAryRef) {
        assert($cntPos==length($accSeqRef->{$acc}));
        my $tmpStr1 = sprintf("%s\t%s\t0\t0\t2\t0\t", ++$ttt, $acc);
        my ($strC1) = ($accSeqRef->{$acc});
        while(my ($chrBefore, $chrAfter) = each(%$eHash1)) {
            $strC1 =~ s|\Q$chrBefore\E|\Q$chrAfter\E \Q$chrAfter\E  |g;
        }
        print $fOutH1 $tmpStr1.$strC1."\n"; #haploview
    }
    for my $i (0..($cntPos-1)) {
        my $loc = $posAryRef->[$i];
        print $fOutH2 join("\t", sprintf("L%04d", $i+1), $loc)."\n"; #haploview_loc
    }
}
sub vntWrite_phase {
    my $self = shift;
    my ($fOut) = @_;
    my ($posAryRef, $accSeqRef) = ($self->{posPickedRef}, $self->{accSeqRef});
    my $fOutH = new IO::File $fOut.".txt", "w";
    my $cntPos = scalar(@$posAryRef);
    my $accAryRef = $self->{accAry};
    my $cntAcc = scalar(@$accAryRef);
    my $eHash2 = {"A"=>"A", "C"=>"C", "G"=>"G", "T"=>"T", "N"=>"?"};
    print $fOutH join("\n", $cntAcc, $cntPos, join("", ("S"x$cntPos)), join(" ", @$posAryRef))."\n"; #phase
    for my $acc (@$accAryRef) {
        assert($cntPos==length($accSeqRef->{$acc}));
        my ($strC2) = ($accSeqRef->{$acc});
        $strC2 =~ s/N/\?/g;
        print $fOutH join("\n", $acc, $strC2, $strC2)."\n"; #phase (cont.)
    }
}
sub vntWrite_rsq {
    my $self = shift;
    my ($fOut) = @_;
    my ($posAryRef, $accSeqRef) = ($self->{posPickedRef}, $self->{accSeqRef});
    my $fOutH = new IO::File $fOut.".txt", "w";
    my $cntPos = scalar(@$posAryRef);
    my $accAryRef = $self->{accAry};
    my $cntAcc = scalar(@$accAryRef);
    my $ttt = 0;
    for my $acc (@$accAryRef) {
        assert($cntPos==length($accSeqRef->{$acc}));
        print $fOutH join("\n", ">".$acc, $accSeqRef->{$acc})."\n\n"; #rsq
    }
}
sub vntWrite_LDhat {
    my $self = shift;
    my ($fOut) = @_;
    my ($posAryRef, $accSeqRef) = ($self->{posPickedRef}, $self->{accSeqRef});
    my $fOutH1 = new IO::File $fOut.".txt", "w";
    my $fOutH2 = new IO::File $fOut."_loc.txt", "w";
    my $cntPos = scalar(@$posAryRef);
    my $accAryRef = $self->{accAry};
    my $cntAcc = scalar(@$accAryRef);
    my ($locMin, $locMax) = (min(@$posAryRef), max(@$posAryRef));
    assert($locMin==$posAryRef->[0]);
    print $fOutH2 sprintf("%d\t%d\tL", $cntPos, ceil(($locMax-$locMin+1)/1000))."\n"; #LDhat_loc
    for my $i (0..($cntPos-1)) {
        my $loc = $posAryRef->[$i];
        print $fOutH2 sprintf("%.3f", ($loc-$locMin+1)/1000)."\n"; #LDhat_loc (cont.)
    }
    print $fOutH1 sprintf(" %d %d 1", $cntAcc, $cntPos)."\n"; #LDhat
    for my $acc (@$accAryRef) {
        assert($cntPos==length($accSeqRef->{$acc}));
        my ($strC2) = ($accSeqRef->{$acc});
        $strC2 =~ s/N/\?/g;
        print $fOutH1 join("\n", ">".$acc, $strC2)."\n"; #LDhat (cont.)
    }
}
sub vntWrite_LDhot {
    my $self = shift;
    my ($fOut) = @_;
    my ($posAryRef, $accSeqRef) = ($self->{posPickedRef}, $self->{accSeqRef});
    my $fOutH1 = new IO::File $fOut.".txt", "w";
    my $fOutH2 = new IO::File $fOut."_loc.txt", "w";
    my $cntPos = scalar(@$posAryRef);
    my $accAryRef = $self->{accAry};
    my $cntAcc = scalar(@$accAryRef);
    my ($locMin, $locMax) = (min(@$posAryRef), max(@$posAryRef));
    assert($locMin==$posAryRef->[0]);
    print $fOutH2 sprintf("%d\t%d\tL", $cntPos, ceil(($locMax-$locMin+1)/1000))."\n"; #LDhat_loc
    for my $i (0..($cntPos-1)) {
        my $loc = $posAryRef->[$i];
        print $fOutH2 sprintf("%.3f", ($loc-$locMin+1)/1000)."\n"; #LDhat_loc (cont.)
    }
    print $fOutH1 sprintf(" %d %d 1", $cntAcc, $cntPos)."\n"; #LDhat
    for my $acc (@$accAryRef) {
        assert($cntPos==length($accSeqRef->{$acc}));
        my ($strC2) = ($accSeqRef->{$acc});
        $strC2 =~ s/N/\?/g;
        print $fOutH1 join("\n", ">".$acc, $strC2)."\n"; #LDhat (cont.)
    }
}
1;
__END__
