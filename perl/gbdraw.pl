#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use File::Basename;
use InitPath;
use Common;
use Bio::Graphics;
use Bio::SeqFeature::Generic;
use Data::Dumper;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;

my @pKs = qw/key glyph height bgcolor fgcolor connector desc/;
my $pVs = { 
  mt_30 => [qw/IMGAG_Gene processed_transcript 6 skyblue slateblue solid 1/],
  mt_35 => [qw/IMGAG_Gene processed_transcript 6 skyblue slateblue solid 1/],
  mt_gi => [qw/MtGI-10.0 transcript 6 pink darkgreen dashed 1/],
  mt_35_gi => [qw/MtGI-10.0 transcript 6 pink darkgreen dashed 1/],
  mt_gea => [qw/Mt51k_Array transcript 6 orangered navy dashed 1/],
  mt_35_gea => [qw/Mt51k_Array transcript 6 orangered navy dashed 1/],
  mt_35_defl => [qw/DEFL_CDS transcript 6 orchid slateblue dashed 1/],
  mt_35_crp_hit =>[qw/CRP_hit_support transcript2 6 tbs springgreen solid 0/],
  mt_35_crp_jcvi => [qw/CRPs_built_by_JCVI processed_transcript 6 gold slateblue solid 1/],
  mt_35_crp_model => [qw/CRP_true_models processed_transcript 6 springgreen slateblue solid 1/],
};


sub plotFeaturesByLoc_from_db {
  my ($loc, $ps, $fo, $hilite, $xtracks) = rearrange([qw/loc ps out hilite xtracks/], @_);
  my ($chr, $s, $e) = map {$loc->$_} qw/seq_id start end/;
  my $panel = Bio::Graphics::Panel->new(
    -start=>$s, -end=>$e, -width=>800, 
    -pad_left=>10, -pad_right=>40, -key_style=>'between', -grid=>1);
  my $full_length = Bio::SeqFeature::Generic->new(
    -display_name=>$chr, -start=>$s, -end=>$e);
  $panel->add_track($full_length, -glyph=>'arrow', -bump=>0, -double=>1, -tick=>1, -label=>1);
  for my $p (@$ps) {
    my ($db, $opt, $ld, $conf, $fes) = map {$p->{$_}} qw/db opt ld conf fes/;
    $ld ||= Localdb->new(-db=>$db, -opt=>$opt);
    $fes ||= [ $ld->getFeatures(-loc=>$loc, -types=>$ld->{type}) ];
#    $fes = [grep {$_->source_tag eq "picked"} @$fes] if $db eq "mt_35_crp_hit";
    $conf ||= getTrackConf($db, $p->{opt});
    my ($glyph, $key, $height, $conn, $bg, $fg, $label, $desc) = 
      map {$conf->{$_}} qw/glyph key height connector bgcolor fgcolor label desc/;
    $panel->add_track(
      $fes, -glyph=>$glyph, -key=>$key, -height=>$height, -connector=>$conn,
      -bgcolor=>$bg, -fgcolor=>$fg, -label=>$label, -description=>$desc, 
      -hilite=>sub {return shift->id eq $hilite ? 'yellow' : undef}
    );
  }
  for my $xtrack (@$xtracks) {
    $panel->add_track(
      $xtrack->{fe},
      -glyph      => $xtrack->{glyph},
      -graph_type => $xtrack->{graph_type},
      -point_symbol => $xtrack->{point_symbol},
      -point_radius => $xtrack->{point_radius},
      -scale        => $xtrack->{scale},
      -height    => $xtrack->{height},
      -key       => $xtrack->{key},
      -connector => $xtrack->{connector},
      -bgcolor   => $xtrack->{bgcolor},
      -fgcolor   => $xtrack->{fgcolor},
      -label     => $xtrack->{label},
      -description => $xtrack->{description},
      -min_score => $xtrack->{min_score},
      -max_score => $xtrack->{max_score}
    );
  }
  open(FH, ">$fo") or die "cannot open $fo for write\n";
  print FH $panel->png;
}
sub drawFe {
    my ($fe, $fo, $ps, $padding, $xtracks) =
        rearrange(["fe", "out", "ps", "padding", "xtracks"], @_);
    $padding ||= 0.25;
    my $loc1 = Bio::Location::Simple->new(-seq_id=>$fe->seq_id, -start=>$fe->start, -end=>$fe->end);
    my @locLimits;
    for my $p (@$ps) {
        $p->{ld} ||= Localdb->new(-db=>$p->{db}, -opt=>$p->{opt});
        my $ld = $p->{ld};
        my @fes = $ld->getFeatures(-loc=>$loc1, -types=>$ld->{type});
        $p->{fes} = \@fes;
        push @locLimits, map {$_->start, $_->end} @fes;
    }
    my ($pS, $pE) = ( min(@locLimits), max(@locLimits) );
    my $pLen = $pE - $pS + 1;
    my $s = int($pS - $pLen * $padding);
    my $e = int($pE + $pLen * $padding);
    my $loc2 = Bio::Location::Simple->new(-seq_id=>$fe->seq_id, -start=>$s, -end=>$e);
    drawRange(-loc=>$loc2, -ps=>$ps, -out=>$fo, -hilite=>$fe->id, -xtracks=>$xtracks);
}
sub drawFes {
    my ($dirO, $ref) = rearrange(['dir', 'param'], @_);
    my ($tag, $ps, $optConf) = map {$ref->{$_}} qw/tag opts optconf/;
    $optConf ||= 1;
    $dirO ||= dir($DIR_Misc2, "gbrowse_figs", $tag);
    system("mkdir -p $dirO") unless -d $dirO;
    system("rm $dirO/*");
    my $p = $ps->[-1];
    my $ld = Localdb->new(-db=>$p->{db}, -opt=>$p->{opt});
    my @fes = $ld->getFeatures(-types=>$ld->{type});
    if($optConf == 2) {
        @fes = grep {$_->source_tag eq "picked"} @fes;
    }
    print "preparing to draw ".@fes." features:\n";
    for my $p (@$ps) {
        my ($db, $opt) = map {$p->{$_}} qw/db opt/;
        my $ld = Localdb->new(-db=>$db, -opt=>$opt);
        $p->{ld} = $ld;
    }
    for my $i (0..$#fes) {
        my $fe = $fes[$i];
        my $fImg = file($dirO, $fe->id.".png");
        drawFe(-fe=>$fe, -out=>$fImg, -ps=>$ps);
        printf "\t%4d done\n", $i+1 if ($i+1) % 100 == 0;
    }
}
sub drawCov {
    my ($chr, $start, $end, $width, $db, $fOut, $xtracks, $regexHash) =
        rearrange(["chr", "start", "end", 'width', "db", "out", "xtracks", "hilite"], @_);
    $width ||= 800;
    my $chrName = $chr=~/^\d$/ ? "MtChr$chr" : $chr;
    my $fOutH = new IO::File $fOut, "w";
    my $locStr = sprintf("%s:%s..%s", $chrName, $start, $end);
    my @genes = getFeatures(-loc=>$locStr, -type=>['mRNA']);
    my $panel = Bio::Graphics::Panel->new(
        -start      => $start,
        -end        => $end,
        -width      => $width,
        -pad_left   => 10,
        -pad_right  => 10,
                -key_style  => 'between',
        -grid       => 1
    );
    my $full_length = Bio::SeqFeature::Generic->new(
        -display_name => $chrName,
        -start    => $start,
        -end      => $end
    );
    $panel->add_track(
        $full_length,
        -glyph  => 'arrow',
        -bump   => 0,
        -double => 1,
        -tick   => 1,
        -label  => 1
    );
    $panel->add_track(
        \@genes,
        -glyph       => 'processed_transcript',
        -height      => 6,
        -key         => 'IMGAG_gene',
        -bgcolor     => 'skyblue',
        -fgcolor     => 'slateblue',
        -connector   => 'solid',
        -label       => 0,
        -description => 0,
        -hilite      => sub {
                my ($feTmp) = @_;
                my $note = join(" |", $feTmp->get_tag_values("Note"));
                my $color = undef;
                for my $regex (keys(%$regexHash)) {
                    $color = $regexHash->{$regex} if $note =~ /$regex/i;
                }
                return $color;
            }
    );
    for my $dbOpt (sort keys(%$db)) {
        my @feAry = getFeatures(
            -loc   => $locStr,
            -type  => $db->{$dbOpt}->{type},
            -db    => $db->{$dbOpt}->{biodb}
        );
        my ($trackKey, $bgColor, $fgColor);
        ($trackKey, $bgColor, $fgColor) = ("MTGI-9.0", "pink", "darkgreen") if $dbOpt eq "mt_gi";
        $panel->add_track(
            \@feAry,
            -glyph       => $dbOpt eq "mt_defl_model" ? 'processed_transcript' : 'transcript',
            -height      => 6,
            -key         => $trackKey,
            -connector   => 'dashed',
            -bgcolor     => $bgColor,
            -fgcolor     => $fgColor,
            -label       => 1,
            -description => 1
        );
    }
    for my $xtrack (@$xtracks) {
        $panel->add_track(
            $xtrack->{fe},
            -glyph      => $xtrack->{glyph},
                        -graph_type => $xtrack->{graph_type},
            -point_symbol => $xtrack->{point_symbol},
                        -point_radius => $xtrack->{point_radius},
                        -scale        => $xtrack->{scale},
            -height    => $xtrack->{height},
            -key       => $xtrack->{key},
            -connector => $xtrack->{connector},
            -bgcolor   => $xtrack->{bgcolor},
            -fgcolor   => $xtrack->{fgcolor},
            -label     => $xtrack->{label},
            -description => $xtrack->{description},
            -min_score => $xtrack->{min_score},
            -max_score => $xtrack->{max_score}
        );
    }
    print $fOutH $panel->png;
}
sub makeCovTrack {
    my ($fe, $accAry, $covOpt, $padding, $flagRatio, $refAcc, $opt) =
        rearrange(['fe', 'acc', 'covopt', 'padding', 'ratio', 'ref', 'opt'], @_);
    $opt ||= "log" if $flagRatio == 1;
    $padding ||= 0.25;
    my $covOptHash = {1=>'cov', 2=>'uniqCov'};
    my $fLen = $fe->end - $fe->start + 1;
    my $start = int($fe->start - $fLen * $padding);
    my $end   = int($fe->end   + $fLen * $padding);
    #print join("\t", $start, $end)."\n";
    my $covHash = getCov(-acc=>$accAry, -chr=>$fe->seq_id, -start=>$start, -end=>$end,
        -covopt=>$covOpt, -ratio=>$flagRatio, -ref=>$refAcc);
    my $rst = [];
    my ($min_score, $max_score);
    for my $acc (sort(keys(%$covHash))) {
        my @covArySorted = sort {$a<=>$b} @{$covHash->{$acc}};
        my $idx;
        for my $i (0..$#covArySorted) {
            $idx = $i;
            if($covArySorted[$i] ne "") {
                last if $flagRatio != 1;
                last if $flagRatio == 1 && $covArySorted[$i] != 0;
            }
        }
        #print join("=", $covArySorted[$idx], $covArySorted[$#covArySorted])."\n";
        $min_score ||= $covArySorted[$idx];
        $max_score ||= $covArySorted[$#covArySorted];
        $min_score = $min_score < $covArySorted[$idx] ? $min_score : $covArySorted[$idx];
        $max_score = $max_score > $covArySorted[$#covArySorted] ? $max_score : $covArySorted[$#covArySorted];
    }
    #print join("-", $min_score, $max_score)."\n";
    $min_score ||= 0.05;
    $max_score ||= 20;
    for my $acc (sort(keys(%$covHash))) {
        my @covAry = @{$covHash->{$acc}};
        my ($values, $min, $max) = transform(-data=>\@covAry, -min=>$min_score, -max=>$max_score, -opt=>$opt);
        my $key = join("_", $acc, $covOptHash->{$covOpt});
        $key .= " [$opt($acc/$refAcc)]" if $flagRatio == 1;
        my $track = { glyph=>"xyplot", graph_type=>"points", #point_symbol=>"disc", point_radius=>4,
            height=>40, key=>$key, fgcolor=>'purple', min_score=>$min, max_score=>$max };
        my $feCov = Bio::SeqFeature::Generic->new(-start=>$start, -end=>$end);
        for my $i (0..@$values-1) {
            my $part = Bio::SeqFeature::Generic->new(
                -start=>$start+$i, -end=>$start+$i, -score=>$values->[$i]);
            $feCov->add_SeqFeature($part);
            #print $covHash->{$acc}->[$i]."=";
        }
        $track->{fe} = $feCov;
        push @$rst, $track;
    }
    return $rst;
}
sub transform {
    my ($data, $min_score, $max_score, $option) = rearrange(['data', 'min', 'max', 'opt'], @_);
    my $optHash = {'log'=>1};
    die("Unsupported transform: $option\n") unless exists $optHash->{$option};
    my ($values, $min, $max) = ([], $min_score, $max_score);
    if($option eq "log") {
        for my $value (@$data) {
            if($value > 0) {
                push (@$values, log($value));
            } elsif($value == 0) {
                push (@$values, log($min_score));
            } else {
                die("ratio[$value] must >= 0\n");
            }
        }
        ($min, $max) = (log($min_score), log($max_score));
    }
    return ($values, $min, $max);
}
sub make_range {
      my $feature = shift;
      my $score        = $feature->score;
      my $range_top    = $score + 5*sqrt($score) + rand(50);
      my $range_bottom = $score - 5*sqrt($score) - rand(50);
      my $quartile_top    = $score + 2*sqrt($score) + rand(50);
      my $quartile_bottom = $score - 2*sqrt($score) - rand(50);
      return [$score,$range_bottom,$range_top,$quartile_bottom,$quartile_top];
}

