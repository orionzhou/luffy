#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use FindBin;
use lib $FindBin::Bin;
use File::Path qw/make_path remove_tree/;
use Data::Dumper;
use List::Util qw/min max sum/;
use Common;

my $org = "Athaliana";
#$org = "Mtruncatula_3.5";
my $dir = "/home/youngn/zhoup/Data/misc2/spada.affy/$org";

my $f01 = "$dir/01_probe.tbl";
my $f02 = "$dir/02_probe.fa";
#get_probe_seq("$dir/00.tbl", $f01, $f02);

my $f13 = "$dir/13.gal";
#probe_to_gene($f13, "$dir/21_crp_gs.gtb", "$dir/23.tbl");
#probe_to_gene($f13, "$dir/31_crp_spada.gtb", "$dir/33.tbl");

sub probe_to_gene {
    my ($fp, $fg, $fo) = @_;
    my $tp = readTable(-in=>$fp, -header=>1);
    my $tg = readTable(-in=>$fg, -header=>1);

    open(FH, ">$fo") or die "cannot write to $fo\n";
    print FH join("\t", qw/id set num probes/)."\n";
    for my $i (0..$tg->nofRow-1) {
        my ($id, $chr, $beg, $end) = map {$tg->elm($i, $_)} qw/id chr beg end/;
        my $tps = $tp->match_pattern("\$_->[6] eq '$chr' && \$_->[7] >= $beg && \$_->[8] <= $end");
        my @probes = $tps->col('qId');
        my ($set, $num) = ('', 0);
        if(@probes > 0) {
            my $h;
            for my $probe (@probes) {
                my ($pset, $cnt) = split(/\./, $probe);
                $h->{$pset} ||= 0;
                $h->{$pset} ++;
            }
            my @sets = sort {$h->{$a} <=> $h->{$b}} keys(%$h);
            $set = $sets[-1];
            $num = $h->{$set};
        }
        print FH join("\t", $id, $set, $num, join(",", @probes))."\n";
    }
    close FH;
}
sub get_probe_seq {
    my ($fi, $fo1, $fo2) = @_;
    my ($fhi, $fho);
    open($fhi, "<$fi") or die "cannot read $fi\n";
    open($fho, ">$fo1") or die "cannot write $fo1\n";
    my $seqHO = Bio::SeqIO->new(-file=>">$fo2", -format=>'fasta');
    
    print $fho join("\t", qw/id set seq/)."\n";
    my $h;
    while(<$fhi>) {
        chomp;
        my ($set, $num1, $num2, $num3, $seq) = split "\t";
        $h->{$set} ||= 0;
        $h->{$set} ++;
        my $id = sprintf "%s.%02d", $set, $h->{$set};
        $seqHO->write_seq( Bio::Seq->new(-id=>$id, -seq=>$seq) );
        print $fho join("\t", $id, $set, $seq)."\n";
    }
    $seqHO->close();
    close $fho;
}


