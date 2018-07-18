package VntEffect;
=head1 NAME
my $myEff = VntEffect->new(-fedb=>$feDb, -refdb=>$refDb);
my @test = qw/chr1:2144758:A:g chr1:2144754:T:a chr5:30092376:C:t chr5:30192504:C:t chr4:18158236:A:g/;
for my $a (@test) {
    my ($vntEffectRef, $vntDescRef, $aaLocRef) = $myEff->getEffect(-vnt=>$a, -context=>'');
    for my $i (0..scalar(@$vntDescRef)-1) {
        my ($vntEffect, $vntDesc, $aaLoc) = ($vntEffectRef->[$i], 
            $vntDescRef->[$i], $aaLocRef->[$i]);
        print join("\t", $vntEffect, $vntDesc, $aaLoc)."\n";
    }
}
=cut
use strict; use Init; use Common; use Localdb; use Readfile;
use Bio::Seq; use Bio::SeqIO; use Data::Dumper;
use Bio::Location::Simple; use Bio::Coordinate::Pair;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Math::Round qw/round/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/vntDesc getEffect/;
@EXPORT_OK = qw//;
my $effectCode = {
    1  => "Non-synonymous",
    2  => "Synonymous",
    3  => "Frameshift",
    4  => "Stop lost",
    5  => "Stop gained",,
    6  => "Start lost",
    7  => "Essential splice site",
    8  => "Splice site",
    9  => "Upstream",
    10 => "5' UTR",
    11 => "Intronic",
    12 => "3' UTR",
    13 => "Downstream",
    14 => "Intergenic"
};
my $effectCodeR = { map {$effectCode->{$_} => $_} keys %$effectCode };
sub new {
    my $self = shift;
    my ($refDb, $feDb) = rearrange(["refdb", "fedb"], @_);
    my $ld = Localdb->new(-db=>$feDb, -opt=>'mRNA');
    return bless {ld => $ld, refdb=>$refDb, fedb=>$feDb},
        ref($self) || $self;
}
sub vntDesc {
    my $self = shift;
    my ($fi, $fo) = rearrange(["in", "out"], @_);
    my $t = readTable(-in=>$fi, -header=>1);
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id chr pos alleles effect desc aaLoc/). "\n";
    for my $i (0..$t->nofRow-1) {
        last if $i < 0;
        my ($id, $chr, $pos, $alleles) = map {$t->elm($i, $_)} qw/id chr pos alleles/;
        my ($ref, @vars) = split(":", $alleles);
        for my $var (@vars) {
            my $vntStr = join(":", $chr, $pos, $ref, $var); 
            my ($refE, $refD, $refL) = $self->getEffect(-vnt=>$vntStr, -context=>$id);
            for my $i (0..scalar(@$refD)-1) {
                my ($vntEffect, $vntDesc, $aaLoc) = ($refE->[$i], $refD->[$i], $refL->[$i]);
                print $fh join("\t", $id, $chr, $pos, "$ref:$var", $vntEffect, $vntDesc, $aaLoc)."\n";
            }
        }
        print "$i\n" if $i % 1_000 == 0 && $i > 0;
    }
}
sub splitCmpStr {
    my ($cmpStr) = @_;
    $cmpStr =~ /^\[([\w\-])\|([\w\-])\]\[(\w)\|(\w)\]$/;
    my ($type, $snpStr);
    if(!$3) {
        print "cannot recognize cmpStr: $cmpStr\n";
        exit -1;
    }
    if($3 eq "R" && $4 eq "R") {
        assert($1 eq $2);
        $type = 0;
    } elsif($1 eq $2) {
        $type = 0;
    } elsif($3 eq "D" || $4 eq "D" || $3 eq "I" || $4 eq "I") {
        $type = -1;
    } elsif($3 eq "N" || $4 eq "N") {
        $type = 0;
    } else {
        unless(($3 eq "S" || $4 eq "S") && $1 ne $2) {
            die "error in splitCmpStr - VntEffect.pm\n"
        }
        ($type, $snpStr) = (1, $1.":".$2);
    }
    return ($type, $snpStr);
}
sub getEffect {
    my $self = shift;
    my ($vntStr, $contextId) = rearrange(['vnt', 'context'], @_);
    my ($seqId, $pos, $ref, $var) = split(":",$vntStr);
    my (@effect, @desc, @aaloc);
    my $typeRef = [qw/CDS five_prime_utr three_prime_utr mRNA gene/];
    my ($flagOut, $posTmp, $mRNA) = (0, int($pos), {});
    my $fe; 
    $fe = $self->{ld}->getFeatureByName($contextId) if $contextId;
    if($fe) {
        $mRNA->{$contextId} = $fe;
    } else {
        my $loc = Bio::Location::Simple->new(-seq_id=>$seqId, -start=>$posTmp, -end=>$posTmp);
        my @fes = $self->{ld}->getFeatures(-loc=>$loc, -types=>['mRNA']);
        for my $fe (@fes) {
            $mRNA->{$fe->id} = $fe if $fe->primary_tag eq "mRNA";
        }
    }
    if(scalar(keys(%$mRNA)) == 0) {
        my $locStr2 = $seqId.":".($posTmp-5000)."..".($posTmp+5000);
        my ($oneEffect, $oneDesc) = 
            $self->effectIntergenic(-range=>$locStr2, -pos=>$posTmp, -opt=>1);
        my $oneAaLoc = "-";
        push (@effect, $oneEffect) && push (@desc, $oneDesc) && push (@aaloc, $oneAaLoc);
    } else {
        if(scalar(keys(%$mRNA)) > 1) {
            print "\t".scalar(keys(%$mRNA))." gene models used for $vntStr\n";
        }
        for my $mRNAId (keys(%$mRNA)) {
            my ($fHash, $rHash, @subFe) = $self->getStructure($mRNA->{$mRNAId});
            my ($flag, $tmpDesc, $tag) = (0, "", 0);
            my ($hitFeature, $hitStart, $hitEnd);
            for my $locStr (keys %$rHash) {
                my @tmpAry = split(/\.\./, $locStr);
                my ($start, $end) = (int($tmpAry[0]), int($tmpAry[1]));
                if ($start<=floor($pos) && $end>=floor($pos)) {
                    $tmpDesc = $rHash->{$locStr};
                    $flag = 1;
                    ($hitStart, $hitEnd) = ($start, $end);
                }
            }
            for my $fe (@subFe) {
                $hitFeature = $fe if $fe->start==$hitStart && $fe->stop==$hitEnd;
            }
            if($flag == 0) {
                my ($oneEffect, $oneDesc) =
                    $self->effectIntergenic(-context=>$mRNA->{$mRNAId}, -pos=>$posTmp, -opt=>2);
                my $oneAaLoc = "-";
                push (@effect, $oneEffect) && push (@desc, $oneDesc) && push (@aaloc, $oneAaLoc);
            } else {
                my ($oneEffect, $oneDesc, $oneAaLoc);
                if ($tmpDesc eq "CDS") {
                    die "no feature hit\n" unless $hitFeature;
                    ($oneEffect, $oneDesc, $oneAaLoc) = $self->effectCds($pos, $ref, $var, @subFe);
                    push (@effect, $oneEffect) && push (@desc, $oneDesc) && push (@aaloc, $oneAaLoc);
                    #&effectCds2($pos, $ref, $var, $hitFeature);
                } else {
                    if ($tmpDesc eq "five_prime_UTR") {
                        $oneEffect = 10;
                    } elsif ($tmpDesc eq "three_prime_UTR") {
                        $oneEffect = 12;
                    } elsif ($tmpDesc eq "intron") {
                        my ($d1, $d2) = (abs($pos-$hitStart), abs($pos-$hitEnd));
                        my $d = $d1<$d2 ? $d1:$d2;
                        if ($d<=1) {
                            $oneEffect = 7;
                        } elsif ($d>=2 and $d<=7) {
                            $oneEffect = 8;
                        } else {
                            $oneEffect = 11;
                        }
                    }
                    $oneDesc = $mRNAId;
                    $oneAaLoc = "-";
                    push (@effect, $oneEffect) && push (@desc, $oneDesc) && push (@aaloc, $oneAaLoc);
                }
                if ($tmpDesc eq "CDS" || $tmpDesc eq "five_prime_UTR" || $tmpDesc eq "three_prime_UTR") {
                    my $aryRef = $self->effectExon($pos, $fHash->{"intron"});
                    if ($aryRef) {
                        ($oneEffect, $oneDesc) = @$aryRef;
                        $oneAaLoc = "-";
                        push (@effect, $oneEffect) && push (@desc, $oneDesc) && push (@aaloc, $oneAaLoc);
                    }
                }
            }
        }
    }
    my @effectMore;
    for my $t (@effect) {
        push (@effectMore, $effectCode->{$t});
    }
    return (\@effectMore, \@desc, \@aaloc);
}
sub effectExon {
    my $self = shift;
    my ($posG, $posIntronLstRef) = @_;
    my $aryRef;
    my $d2splice = 100;
    for my $posIntron (@$posIntronLstRef) {
        my ($startI, $endI) = ($posIntron->[0], $posIntron->[1]);
        my ($d1, $d2) = (abs($startI-$posG), abs($endI-$posG));
        my $d2spliceHere = $d1<$d2 ? $d1:$d2;
        $d2splice = $d2splice<$d2spliceHere ? $d2splice:$d2spliceHere;
    }
    if ($d2splice <= 3) {
        $aryRef = [8, $d2splice." bp into exon"];
    }
    return $aryRef;
}
sub effectIntergenic {
    my $self = shift;
    my ($rangeLocStr, $pos, $context, $option) = rearrange(['range', 'pos', 'context', 'opt'], @_);
    my $optHash = {1=>'by range', 2=>'by feature context'};
    die("unsupported option: $option\n") unless exists $optHash->{$option};
    die("range LocStr must be specified\n") if $option == 1 && !$rangeLocStr;
    die("context must be provided\n") if $option == 2 && !$context;
    $pos = int($pos);
    my ($effect, $desc);
    if($option == 1) {
        my ($sid, $s, $e) = splitLocStr($rangeLocStr);
        my $loc = Bio::Location::Simple->new(-seq_id=>$sid, -start=>$s, -end=>$e);
        my @sFeatures = $self->{ld}->getFeatures(-loc=>$loc, -types=>["mRNA"]);
        if (@sFeatures == 0) {
            $effect = 14;
            $desc = "-";
        } else {
            my $featureN;
            my $distanceN = 5000;
            for my $sFe (@sFeatures) {
                my ($d1, $d2) = (abs($sFe->start-$pos), abs($sFe->end-$pos));
                my $distance = $d1<$d2 ? $d1:$d2;
                if ($distance < $distanceN) {
                    $distanceN = $distance;
                    $featureN = $sFe;
                }
            }
            if (( $featureN->end < $pos && $featureN->strand == 1 ) ||
                    ( $featureN->start > $pos && $featureN->strand == -1)) {
                $effect = 13;
                $desc = "$distanceN bp";
            } else {
                $effect = 9;
                $desc = "$distanceN bp";
            }
        }
    } else {
        my $fe = $context;
        my ($d1, $d2) = (abs($fe->start-$pos), abs($fe->end-$pos));
        my $distance = $d1<$d2 ? $d1:$d2;
        if (( $fe->end < $pos && $fe->strand == 1 ) ||
                ( $fe->start > $pos && $fe->strand == -1)) {
            $effect = 13;
            $desc = "$distance bp";
        } else {
            $effect = 9;
            $desc = "$distance bp";
        }
    }
    return ($effect, $desc);
}
sub effectCds {
    my $self = shift;
    my ($posG, $ref, $var, @subFe) = @_;
    my ($effect, $desc, $aaloc);
    my (%cdsHash, $strand);
    my ($cdsLen, $cds, $posL) = (0, "", 0);
    for my $fe (@subFe) {
        if ($fe->primary_tag eq "CDS") {
            $cdsHash{$fe->start} = $fe;
            $strand = $fe->strand if !$strand;
            die "strand not consistent" unless $strand eq $fe->strand;
        }
    }
    print $subFe[0]->id." => ".$cdsLen." bp => error\n" if $cdsLen % 3 != 0;
    my @cdsStarts = sort {$a<=>$b} keys %cdsHash;
    @cdsStarts = reverse sort {$a<=>$b} keys %cdsHash if $strand == -1;
    for my $cdsStart (sort keys %cdsHash) {
        my $fe = $cdsHash{$cdsStart};
        my $loc = Bio::Location::Simple->new(-seq_id=>$fe->seq_id, -start=>$fe->start, -end=>$fe->end, -strand=>$strand);
        $cds .= seqRet($loc, $self->{refdb});
        if ($posG>=$fe->start && $posG<=$fe->end) {
            $posL = $posG - $fe->start + 1 + $cdsLen;
            $posL = $fe->end - $posG   + 1 + $cdsLen if $strand == -1;
        }
        $cdsLen += $fe->end - $fe->start + 1;
    }
    my ($cntAA, $posL_AA, $offSet) = ($cdsLen/3, int(($posL-1)/3), ($posL-1)%3);
    if($ref eq "-") {
        if (length($var) % 3 != 0) {
            $effect = 3;
        } else {
            $effect = 1;
        }
        $desc = length($ref)."bp ins";
    } elsif($var eq "-") {
        if (length($ref) % 3 != 0) {
            $effect = 3;
        } else {
            $effect = 1;
        }
        $desc = length($ref)."bp del";
    } else {
        if ($strand == -1) {
            my $refSeq = Bio::Seq->new(-seq=>$ref);
            $ref = $refSeq->revcom->seq;
            my $varSeq = Bio::Seq->new(-seq=>$var);
            $var = $varSeq->revcom->seq;
        }
        my $refAbs = substr($cds, $posL_AA*3+$offSet, 1);
        if($var eq $refAbs) {
            my $tmp = $ref;
            ($ref, $var) = ($var, $tmp);
        }
        my $codonRef = substr($cds, $posL_AA*3, 3);
        if($ref ne substr($cds, $posL_AA*3+$offSet, 1)) {
            print "Non-ref substitution: $posG [$refAbs] $ref -> $var\n";
            $codonRef = substr($cds, $posL_AA*3, $offSet).$ref.substr($cds, $posL_AA*3+$offSet+1, 2-$offSet);
        }
        my $codonVar = substr($cds, $posL_AA*3, $offSet).$var.substr($cds, $posL_AA*3+$offSet+1, 2-$offSet);
        my ($aaRef, $aaVar);
        if ($codonRef =~ /[^ATCG]/) {
            ($aaRef, $aaVar) = ("x", "x");
        } else {
            my $codonRefSeq = Bio::Seq->new(-seq=>$codonRef);
            my $codonVarSeq = Bio::Seq->new(-seq=>$codonVar);
            ($aaRef, $aaVar) = ($codonRefSeq->translate->seq, 
                $codonVarSeq->translate->seq);
        }
        if ($aaRef eq $aaVar) {
            $effect = 2;
        } elsif($aaRef eq "*" && $aaVar ne "*") {
            $effect = 4;
        } elsif($aaRef ne "*" && $aaVar eq "*") {
            $effect = 5;
        } elsif($posL_AA == 0) {
            $effect = 6;
        } else {
            $effect = 1;
        }
        my $codonChange = join("->",$codonRef, $codonVar);
        my $aaChange = join("->",$aaRef, $aaVar);
        $desc = "$codonChange [$aaChange] [$posL/$cdsLen]";
    }
    $aaloc = ($posL_AA+1)."/".$cntAA;
    return ($effect, $desc, $aaloc);
}
sub effectCds2 {
    my $self = shift;
    my ($pos, $ref, $var, $fe) = @_;
    my ($effect, $desc, $aaloc);
    my $cds = $fe->seq->seq;
    my ($hitStart, $hitEnd) = ($fe->start, $fe->end);
    if ($fe->strand == -1) {
        my $seqRef = Bio::Seq->new(-seq=>$ref);
        $ref = $seqRef->revcom->seq;
        my $seqVar = Bio::Seq->new(-seq=>$var);
        $var = $seqVar->revcom->seq;
    }
    my $phase = $fe->phase;
    print $phase."\t";
    my $global = Bio::Location::Simple->new
        (-seq_id=>"chr", -start=>$hitStart, -end=>$hitEnd, -strand=>1);
    my $local = Bio::Location::Simple->new
        (-seq_id=>"cds", -start=>1, -end=>$hitEnd-$hitStart+1, -strand=>$fe->strand);
    my $pair = Bio::Coordinate::Pair->new(-in=>$global, -out=>$local);
    my $posSys1 = Bio::Location::Simple->new(-start=>$pos, -end=>$pos);
    my $posSys2 = $pair->map($posSys1);
    my $posL = $posSys2->start;
    print $ref."/".substr($cds, $posL-1, 1)."\t";
    my $posL_AA = int(($posL-$phase-1)/3);
    my $offSet = ($posL-$phase-1) % 3;
    my $refCodon = substr($cds, $phase+3*$posL_AA, 3);
    my $varCodon = substr($cds, $phase+3*$posL_AA, $offSet).$var.
        substr($cds, $phase+3*$posL_AA+$offSet+1, 2-$offSet);
    print join("->", $refCodon, $varCodon)."\n";
    return ($effect, $desc, $aaloc);
}
sub getStructure {
    my $self = shift;
    my ($mRNA) = @_;
    my $rHash = {};
    my ($posAryCds, $posAryMRNA, $posAryUtr5, $posAryUtr3) = 
        ([], [], [], [], []);
    my $posAryExon = [];
    push(@$posAryMRNA, [$mRNA->start, $mRNA->stop]);
    my @subFe = $mRNA->get_SeqFeatures();
    my @tags = qw/exon CDS five_prime_UTR three_prime_UTR intron/;
    for my $fe (@subFe) {
        #print join("\t", $fe->id, $fe->start, $fe->end)."\n";
        my $tag = first_index {$_ eq $fe->primary_tag} @tags;
        if($tag == 0) {
            push(@$posAryExon, [$fe->start, $fe->stop]);
        } elsif($tag == 1) {
            push(@$posAryCds, [$fe->start, $fe->stop]);
        } elsif($tag == 2) {
            push(@$posAryUtr5, [$fe->start, $fe->stop]);
        } elsif($tag == 3) {
            push(@$posAryUtr3, [$fe->start, $fe->stop]);
        }
        if($tag >= 1 && $tag <= 3) {
            $rHash->{$fe->start."..".$fe->stop} = $fe->primary_tag;
        }
    }
    my ($posAryIntron, $lenIntron) = getPosDiff($posAryMRNA, $posAryExon);
    for my $ref (@$posAryIntron) {
        $rHash->{join("..", @$ref)} = "intron";
    }
    my $fHash = {"five_prime_UTR"=>$posAryUtr5, "CDS"=>$posAryCds, 
        "intron"=>$posAryIntron, "three_prime_UTR"=>$posAryUtr3};
    return ($fHash, $rHash, @subFe);
}
1;
__END__
