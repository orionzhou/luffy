#!/usr/bin/perl -w
use strict;
use Common;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use LWP;
use LWP::UserAgent;
use HTML::TableExtract;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
@EXPORT = qw/postUrl fetch_plaza fetch_mtgi/;
@EXPORT_OK = qw//;
sub postUrl {
    my( $url, $formref ) = @_;
    # set up a UserAgent object to handle our request
    my $ua = new LWP::UserAgent(timeout => 300);
    $ua->agent('perlproc/1.0');
    my $response = $ua->post( $url, $formref );
    if( $response->is_success ){
        return $response->content;
    } else {
        return undef;
    }
}
sub fetch_plaza {
    my $DIR_Circos;
    my $url = "http://bioinformatics.psb.ugent.be/legume-plaza/genes/from";
    my $param = {
#  "mt" => ["Medicago%20truncatula", 2330, 46587],
#  "gm" => ["Glycine%20max", 2320, 46382],
#  "lj" => ["Lotus%20japonicus", 1514, 30266],
        "vv" => ["Vitis%20vinifera", 1318, 26346]
        };
    my @hs = qw/Gene Transcript Start Stop Strand Chr Type/;
    push @hs, "Gene family";
    for my $sp (keys %$param) {
        my $te = HTML::TableExtract->new( headers => \@hs );
        my ($spd, $numP, $numR) = @{$param->{$sp}};
        my $dirT = "/tmp/$sp";
#  unlink $dirT;
        my $fo = "$DIR_Circos/02_coords/$sp.txt";
        my $fh = new IO::File $fo, "w";
        print $fh join("\t", qw/gene transcript start end strand chr type family/)."\n";
        my $iStr = join(",", (1..$numP));
        my $cmd = "wget $url/$spd/page:{$iStr} -P $dirT --user=medicago --password=medicago*rules";
#  system($cmd);
        for my $i (1..$numP) {
            my $f = file($dirT, "page\:$i");
            $te->parse_file(new IO::File $f, "r");
        }
        $sp =~ s/^(\w)/\U$1/;
        for my $ts ($te->tables) {
            print $ts->coords."\n";
            for my $row ($ts->rows) {
                my @rows = grep s/(^\s*)|(\s*$)//g, @$row;
                my $chr = $rows[5];
                if($chr =~ /^(\d+)$/ || $chr =~ /^chr(\d+)/i || $chr =~ /\Q$sp\E(\d+)/i) {
                    $rows[5] = "$sp$1";
                    $rows[5] = sprintf("%s%02d", $sp, $1) if $sp eq "Vv";
                }
                print $fh join("\t", @rows)."\n";
            }
        }
    }
}
sub fetch_mtgi {
    my ($fi, $fo) = rearrange(['in', 'out'], @_);
#    $fi = "$DIR_Misc2/mtgi/lib.txt";
#    $fo = "$DIR_Misc2/mtgi/lib_tc.txt";
    my $fhi = new IO::File $fi, "r";
    my $fho = new IO::File $fo, "w";
    my $cnt = 0;
    print $fho join("\t", qw/lib_id lib_name supplier tissue organ disease stage no_in_tc
        no_singleton no_est tc_in_lib tc/), "\n"; 
    my $cntTc = 0;
    while( <$fhi> ) {
        chomp;
        next if !$_ || $_=~/^lib\_id/;
        last if ++$cnt < 1;
        my ($libId, $libName, $supplier, $accession, $tissue, $organ,
            $disease, $stage, $no_in_tc, $no_singleton, $no_est)
            = split("\t", $_);
        $no_in_tc = $no_in_tc ? $no_in_tc : 0;
        $no_singleton = $no_singleton ? $no_singleton : 0;
        my $url1 = "http://compbio.dfci.harvard.edu/cgi-bin/tgi/lib_report.pl?CATNUM=%s&species=Medicago";
        my $url = sprintf($url1, uri_escape($libId));
        my $html = get($url) || die("cannot get html\n");
        my @tcId;
        while($html =~ /\<OPTION\>(\w+)/ig) {
            push (@tcId, $1);
        }
        my $cntTc_in_lib = scalar(@tcId);
        if($html =~ /(\d+) TCs have ESTs/) {
            die "TC number incorrect: $1 != $cntTc_in_lib" unless $cntTc_in_lib == $1;
        } else {
            print "error in parsing html\n";
        }
        print $fho join("\t", $libId, $libName, $supplier, $accession, $tissue, $organ,
            $disease, $stage, $no_in_tc, $no_singleton, $no_est,
            $cntTc_in_lib, join(",", @tcId))."\n";
        print join("\t", $libId, $libName, $cntTc_in_lib),"\n";
        $cntTc += $cntTc_in_lib;
    }
    print "\n", "$cntTc TCs in total\n";
}


1;
__END__
