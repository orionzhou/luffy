#!/usr/bin/perl -w
use strict;
use Bio::Graphics;
use Bio::SeqFeature::Generic;
use Data::Dumper;
use List::Util qw/min max sum/;

use Common;
use BioFeature;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/plot_transcripts_by_loc/;
@EXPORT_OK = qw//;

our @color_configs = (
    [qw/skyblue slateblue/],
    [qw/pink darkgreen/],
    [qw/orangered navy/],
    [qw/orchid slateblue/],
    [qw/orange springgreen/],
    [qw/gold slateblue/],
    [qw/springgreen slateblue/],
);

sub plot_transcripts_by_loc {
    my ($chr, $beg, $end, $fo, $p) = rearrange([qw/chr beg end out p/], @_);
    
    my $panel = Bio::Graphics::Panel->new( -start=>$beg, -end=>$end, -width=>800, -key_style=>'between', -grid=>1 );
    
    my $ruler = Bio::SeqFeature::Generic->new( -display_name=>$chr, -start=>$beg, -end=>$end );
    $panel->add_track($ruler, -glyph=>'arrow', -double=>1, -tick=>2, -label=>1);
 
    for my $i (1..@$p) {
        my ($key, $tg) = @{$p->[$i-1]};
        my $fes = gtb2Features($tg, $chr, $beg, $end);
        
        my ($bg, $fg) = @{ $color_configs[($i % @color_configs) - 1] };
        $panel->add_track($fes, -glyph=>'processed_transcript', 
            -bgcolor=>$bg, -fgcolor=>$fg, -connector=>'solid',
            -implied_utrs=> 1, -label=>1, -description => 0, -key=>$key);
    }
    
    open(FH, ">$fo") or die "cannot open $fo for writing\n";
    print FH $panel->png;
    close FH;
}

1;
__END__

