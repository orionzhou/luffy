package Annotate;
use strict; use Init; use Common; use Localdb; use Mapping; 
use Path::Class; use DBI; 
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/getAnn annForNBS annProbe annHclust annFromGb/;
@EXPORT_OK = qw//;
sub getAnn {
    my ($fIn, $header) = (rearrange['in', 'head'], @_);
    die("$fIn is not there\n") unless -s $fIn;
    $header ||= 0;
    my $fInH = new IO::File $fIn, "r";
    my $firstLine = readline($fInH) if $header;
    my $annH = {};
    while( <$fInH> ) {
        chomp;
        next unless $_;
        my @eleAry = split("\t", $_);
        $annH->{$eleAry[0]} = [] unless exists $annH->{$eleAry[0]};
        push @{$annH->{$eleAry[0]}}, \@eleAry;
    }
    return $annH;
}
sub annForNBS {
    my ($fi, $fo, $db, $idAry, $opt) = rearrange(['in', 'out', 'db', 'ids', 'opt'], @_);
    my $annH = getAnn(-in=>$fi, -head=>1);
    $opt ||= 1;
    my $fh = new IO::File $fo, "w";
    my $ld = Localdb->new(-db=>$db);
    for my $id (@$idAry) {
        my $fe = $ld->getFeatureByName($id);
        die "no feature called $id\n" unless $fe;
        my $loc = Bio::Location::Simple->new(-seq_id=>$fe->seq_id, -start=>$fe->start, -end=>$fe->end);
        my @fes2 = $ld->getFeatures(-loc=>$loc, -types=>'mRNA');
        my @anns;
        for my $fe2 (@fes2) {
            my $id2 = $fe2->id;
            $id2 =~ s/\_[a-zA-Z0-9]+$//;
            my $ann = $id2;
            for my $aryRef (@{$annH->{$id2}}) {
                $ann .= "[".$aryRef->[3]."]";
            } 
            push @anns, $ann;
        }
        my $annStr = join(";", @anns);
        if($opt == 1) {
            print $fh join("\t", $id, $annStr)."\n";
        } elsif($opt == 2) {
            print $fh join("\t", $id, "id {$id}", "clade {$annStr}")."\n";
        }
    }
}
sub annProbe {
    my ($fOut, $idAry, $opt) = rearrange(['out', 'ids', 'opt'], @_);
    $opt ||= 1;
    my $fOutH = new IO::File $fOut, "w";
    my $ld = Localdb->new(-db=>'mt_gea');
    for my $id (@$idAry) {
        my $fe = $ld->getFeatureByName($id);
        die "no feature called $id\n" unless $fe;
        my $loc = Bio::Location::Simple->new(-seq_id=>$fe->seq_id, -start=>$fe->start, -end=>$fe->end);
        my @fes2 = $ld->getFeatures(-loc=>$loc, -types=>'mRNA');
        my @anns;
        for my $fe2 (@fes2) {
            my $id2 = $fe2->id;
            $id2 =~ s/\_[A-Z]$//;
            push @anns, $id2;
        }
        if($opt == 1) {
            print $fOutH join("", map {join("\t", $id, $_)."\n"} @anns);
        } elsif($opt == 2) {
            my $annStr = join(";", @anns);
            print $fOutH join("\t", $id, "probe {$annStr}")."\n";
        }
=exp
        my @anns;
        my @feIdAry = map {$_->id} @fes;
        my ($exp, $condS, $condD) = getExpression(@feIdAry);
        my @idHitAry = sort(keys(%$exp));
        for my $i (0..@idHitAry-1) {
            my $vName = "probeSet".($i+1);
            my $expValue = $exp->{$idHitAry[$i]};
            my @posX = (1..@$expValue);
            my $scaledX = scaleNumber(-value=>\@posX, -down=>0, -up=>0.2);
            my $scaledY = scaleNumber(-value=>$expValue, -down=>0, -up=>0.2);
            my @v = pairwise {($a, $b)} @$scaledX, @$scaledY;
            push @anns, join(" ", $vName, "{".join(" ", @$scaledY)."}");
        }
        print $fOutH join("\t", $id, join(";", @anns))."\n";
=cut
    }
}
sub annHclust {
    my ($fIn1, $fIn2, $fOut) = rearrange(['in1', 'in2', 'out'], @_);
    my $annH1 = getAnn(-in=>$fIn1, -head=>1);
    my $annH2 = getAnn(-in=>$fIn2, -head=>1);
    my $fOutH = new IO::File $fOut, "w";
    for my $geneid (keys %$annH1) {
        my $aryRef1 = $annH1->{$geneid};
        my @probeids = map {$_->[1]} @$aryRef1;
        my @anns;
        for my $probeid (@probeids) {
            next unless exists $annH2->{$probeid};
            my $aryRef2 = $annH2->{$probeid};
            my @annTmp = map {$_->[1]} @$aryRef2;
            push @anns, @annTmp;
        }
        my $annStr = join(";", @anns);
        print $fOutH join("\t", $geneid, "hclust {$annStr}")."\n";
    }
}
sub annFromGb {
    my ($fIn, $fOut1, $fOut2) = @_;
    my $seqInH = Bio::SeqIO->new(-file => $fIn, -format => 'genbank');
    my $seqOutH = Bio::SeqIO->new(-file => ">".$fOut1, -format => 'fasta');
    my $fOutH2 = $fOut2->open("w") or die "Can't create $fOut2: $!";
    my (@lst_orgn, @lst_acc, @lst_gi, @lst_id, %count, %map);
    while(my $seq = $seqInH->next_seq()) {
        for my $feat_obj ($seq->get_SeqFeatures) {        
            if ($feat_obj->has_tag("organism")) {
                my @tmp = $feat_obj->get_tag_values("organism");
                my $orgn = $tmp[0];
                $orgn =~ s/^\W*(.*)\W*$/$1/g;
                my @wordAry = split(/[\W]/,$orgn);
                $orgn = join("_", @wordAry);
                $count{$orgn} = 0 unless exists $count{$orgn};
                $count{$orgn} ++;
                my $id = join("", map {substr($_, 0, 1)} @wordAry);
                my $i = 0;
                while(exists $map{$id}) {
                    $id .= substr($wordAry[-1], ++$i, 1);
                }
                $map{$id} = $orgn;
                $id .= "_".sprintf("%02d",$count{$orgn});
                push @lst_orgn, $orgn;
                push @lst_id, $id;
                my $newseq = Bio::Seq->new(-id=>$id, -seq=>$seq->seq());
                $seqOutH->write_seq($newseq);
            }
        }
        push @lst_acc, $seq->accession_number;
        push @lst_gi, $seq->primary_id;
    }
    for(my $i=0; $i<@lst_orgn; $i++) {
        print $fOutH2 join("\t", $lst_id[$i],"acc {$lst_acc[$i]}",
            "gi {$lst_gi[$i]}", "orgn {$lst_orgn[$i]}", "\n");
    }
    foreach my $orgn (reverse sort {$count{$a}<=>$count{$b}} keys(%count)) {
        print join("\t", $count{$orgn}, $orgn)."\n";
    }
}


1;
__END__
