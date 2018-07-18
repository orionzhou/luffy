package CircosConf;
use strict; use Init; use Common; use Path::Class; use IO::File; use File::Copy;
use Data::Dumper; use Text::Template;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use LWP::UserAgent;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw//;
@EXPORT_OK = qw/write setChr addZoom addPlot addLink setLinks/;
=usage
my $cc = CircosConf->new(-sp=>'mt');
$cc->setChr(-optsp=>2);
$cc->addZoom(-chr=>"chr5", -start=>'8u', -end=>'9u', -scale=>5);
$cc->addZoom(-chr=>"chr9", -start=>'9.1u', -end=>'10.4u', -scale=>5);
$cc->addPlot(-type=>'histogram', -file=>'jkj.txt', -r0=>'0.8r', -r1=>'0.9r', -min=>0, -max=>1, -color_fg=>"red", -color_bg=>"vlgrey");
$cc->setLinks(-radius=>'0.8r');
$cc->addLink(-file=>"link01.txt", -ribbon=>1, -rules=>[[1, 'grey', 1], ["_CHR1_ eq 'chr5'", 'red', 2]]);
$cc->write();
=cut
sub new {
    my $self = shift;
    my ($fo) = rearrange(['out'], @_);
    my $dirW = $DIR_Circos;
    my $fSrc = file($dirW, "01_conf", "00.tmpl");
    my $tmpl = Text::Template->new(-type=>'FILE', -source=>$fSrc);
    my $hash = {fo=>$fo};
    return bless {tmpl=>$tmpl, hash=>$hash, dirW=>$dirW},
        ref($self) || $self;
}
sub setChr {
    my $self = shift;
    my ($optSp, $chrShow, $chrs) = rearrange(['optsp', 'chrshow', 'chrs'], @_);
    $optSp ||= 1;
    $chrShow ||= 0;
    $self->{hash}->{chr_show} = $chrShow;
    my $dirC = dir($self->{dirW}, "01_conf");
    my $hSp = {1=>'rosids', 2=>'assembly_mt_35', 3=>'assembly_mt_30'};
    my $sp = $hSp->{$optSp};
    die "unknown species opt: $optSp\n" unless $sp;
    $self->{hash}->{file_karyotype} = "$dirC/$sp.txt";
    if($chrs) {
        $self->{hash}->{chrs} = join(";", @$chrs);
        $self->{hash}->{chrs_order} = join(";", @$chrs);
    } else {
        if($optSp == 2 || $optSp == 3) {
            $self->{hash}->{chrs} = join(";", map {"chr$_"} (1..8));
            $self->{hash}->{chrs_order} = join(";", map {"chr$_"} (1..8));
        }
    }
}
sub addPlot {
    my $self = shift;
    my ($type, $file, $r0, $r1, $min, $max, $color_fg, $color_bg) = rearrange([qw/type file r0 r1 min max color_fg color_bg/], @_);
    my @strs = (
        "\t<plot>",
        "\t\ttype = $type",
        "\t\tfile = $file",
        "\t\tr0 = $r0",
        "\t\tr1 = $r1",
        "\t\tmin = $min",
        "\t\tmax = $max",
        "\t\tfill_color = $color_fg",
        "\t\tbackground_color = $color_bg",
        "\t\taxis = yes",
        "\t\taxis_color = lgrey",
        "\t</plot>"
    );
    $self->{plots} ||= [];
    push @{$self->{plots}}, join("\n", @strs);
}
sub setLinks {
    my $self = shift;
    my ($radius, $crest, $radius_b1, $radius_b2, $thickness, $color) = rearrange([qw/radius crest radius_b1 radius_b2 thickness color/], @_);
    $radius ||= "0.75r";
    $crest ||= 0.5;
    $thickness ||= 2;
    $color ||= 'black';
    $radius_b1 ||= "0.2r";
    $radius_b2 ||= 0.5;
    my @strs = (
        "\tradius = $radius", 
        "\tcrest = $crest", 
        "\tthickness = $thickness", 
        "\tcolor = $color", 
        "\tbezier_radius = $radius_b1", 
        "\tbezier_radius_purity = $radius_b2"
    );
    $self->{linkProp} ||= join("\n", @strs);
}
sub addLink {
    my $self = shift;
    my ($tag, $file, $ribbon, $rules) = rearrange(['tag', 'file', 'ribbon', 'rules'], @_);
    my @strs = (
        "\t<link $tag>", 
        "\t\tfile = $file",
        "\t\tribbon = $ribbon",
        "\t\t<rules>",
        "\t\t</rules>",
        "\t</link>"
    );
    my @ruleStrs;
    for my $rule (@$rules) {
        my ($condition, $color, $z) = @$rule;
        my @ruleStrsOne = (
            "\t\t\t<rule>",
            "\t\t\t\tcondition = $condition",
            "\t\t\t\tcolor = $color",
            "\t\t\t\tz = $z",
            "\t\t\t</rule>"
        );
        push @ruleStrs, join("\n", @ruleStrsOne);
    }
    splice @strs, 4, 0, @ruleStrs;
    $self->{links} ||= [];
    push @{$self->{links}}, join("\n", @strs);
}
sub addZoom {
    my $self = shift;
    my ($c, $s, $e, $scale) = rearrange(['chr', 'start', 'end', 'scale'], @_);
    my @strs = (
        "\t<zoom>",
        "\t\tchr = $c",
        "\t\tstart = $s",
        "\t\tend = $e",
        "\t\tscale = $scale",
        "\t</zoom>"
    );
    $self->{zooms} ||= [];
    push @{$self->{zooms}}, join("\n", @strs);
}
sub cleanup {
    my $self = shift;
    my $dirC = dir($self->{dirW}, "01_conf");
    $self->{hash}->{file_config_colors}   ||= "$dirC/colors.conf";
    $self->{hash}->{file_config_fonts}    ||= "$dirC/fonts.conf";
    $self->{hash}->{file_config_ideogram} ||= "$dirC/ideogram.conf";
    $self->{hash}->{file_config_ticks}    ||= "$dirC/ticks_genome.conf";
    $self->setChr(-optsp=>1) unless exists $self->{hash}->{file_karyotype};
    $self->{hash}->{zooms} = join("\n", @{$self->{zooms}}) if exists $self->{zooms};
    $self->{hash}->{plots} = join("\n", @{$self->{plots}}) if exists $self->{plots};
    if(exists $self->{links}) {
        $self->setLinks() unless exists $self->{linkProp};
        $self->{hash}->{links} = join("\n", $self->{linkProp}, @{$self->{links}});
    }
    $self->{hash}->{fo} ||= $self->{dirW}."/tmp.png";
}
sub write {
    my $self = shift;
    my ($fo) = @_;
    $fo ||= "/tmp/tmp.conf";
    $self->cleanup();
    my $tmpl = $self->{tmpl};
    my $tmpl_filled = $tmpl->fill_in(hash=>$self->{hash});
    my $fh = new IO::File $fo, "w";
    print $fh $tmpl_filled;
    close $fh;
    print "circos configuration written to $fo\n";
}

1;
__END__
