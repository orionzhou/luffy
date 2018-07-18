#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use InitPath;
use Common;
use Bio::Seq; use Bio::SeqIO; use Graph; use Bio::SeqFeature::Generic;
use Readfile; use Writefile; use Annotate; use Align;
use Time::HiRes qw/gettimeofday tv_interval/; use Data::Dumper;use Path::Class; 
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Math::Round qw/round/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
my $ps = { 
      1 => {tag=>'01', ks=>[2..10], nl=>412, ni=>474, np=>34, reps=>[1..5]},
      2 => {tag=>'02', ks=>[2..10], nl=>300, ni=>443, np=>32, reps=>[1..5]},
      3 => {tag=>'03', ks=>[2..10], nl=>412, ni=>474, np=>34, reps=>[1]},
    11 => {tag=>'11_acc289', ks=>[2..10], nl=>5000, ni=>289, np=>289, reps=>[1]},
    12 => {tag=>'12_acc261', ks=>[2..10], nl=>5000, ni=>261, np=>261, reps=>[1]},
};
my $p = $ps->{12};
#run_structure($p);
#sum_structure($p);
#plot_structure($p);
sub run_structure {
    my ($p) = @_;
    my ($tag, $nl, $np, $ni) = map {$p->{$_}} qw/tag nl np ni/;
    my $dirW = dir($DIR_Misc3, "structure", $tag);
    chdir($dirW) or die "cannot cd to $dirW\n";
    my $f01 = file($dirW, "01_in.txt");
    my $d02 = dir( $dirW, "02_structure");
    system("mkdir -p $d02") unless -d $d02;
    for my $k (@{$p->{ks}}) {
        for my $i (@{$p->{reps}}) { 
            my $f02 = file($d02, "k_$k\_rep_$i");
            runCmd("structure -i $f01 -L $nl -N $ni -K $k -o $f02", 1);
        }
    }
}
sub sum_structure {
    my ($p) = @_;
    my ($tag, $nl, $np, $ni) = map {$p->{$_}} qw/tag nl np ni/;
    my $dirW = dir($DIR_Misc3, "structure", $tag);
    my $f01 = file($dirW, "01_in.txt");
    my $d02 = dir( $dirW, "02_structure");
    my $f03 = file($dirW, "03_prob.txt");
    my $d04 = file($dirW, "04_qmatrix");
    my $d05 = file($dirW, "05_stat");
    system("mkdir -p $d04") unless -d $d04;
    system("mkdir -p $d05") unless -d $d05;
    my $fh3 = new IO::File $f03, "w";
    print $fh3 join("\t", qw/K Run P(D) Mean[P(D)] Var[P(D)] Mean[alpha]/)."\n";
    for my $k (@{$p->{ks}}) {
        for my $i (@{$p->{reps}}) { 
            my $f02 = file($d02, "k_$k\_rep_$i\_f");
            next unless -s $f02;
            my ($prob, $q1, $q2) = parse_structure($f02);
            print $fh3 join("\t", $k, $i, @$prob)."\n";
            die "not $np pops: ".@$q1." in $f02\n" unless $np == @$q1;
            die "not $ni inds: ".@$q2." in $f02\n" unless $ni == @$q2;
            my $pre = "k_$k\_rep_$i";

            my $f041 = file($d04, "$pre\_pop.txt");
            my $fh41 = new IO::File $f041, "w";
            print $fh41 join("\n", @$q1)."\n";
            my $f042 = file($d04, "$pre\_ind.txt");
            my $fh42 = new IO::File $f042, "w";
            print $fh42 join("\n", @$q2)."\n";
            
            my $f05 = file($d05, "$pre.tbl");
            my $fh5 = new IO::File $f05, "w";
            print $fh5 join("\t", qw/id pop k mean up95 dw95/)."\n";
            for my $line (@$q2) {
                my @ps = split(" ", $line);
                my ($id, $pop) = @ps[1,3];
                for my $x (1..$k) {
                    my $mean = $ps[4+$x];
                    my ($up, $dw);
                    if( $ps[4+$k+$x] =~ /^\(([\d\.]+)\,([\d\.]+)\)$/) {
                        ($up, $dw) = ($1, $2);
                    } else {
                        die "unknown element: $ps[$i]\n";
                    }
                    print $fh5 join("\t", $id, $pop, $x, $mean, $up, $dw)."\n";
                }
            }
            print join("\t", "K=$k", "Rep=$i", "Okay")."\n";
        }
    }
}
sub parse_structure {
    my ($fi) = @_;
    my $fh = new IO::File $fi, "r";
    my $prob = [('') x 6];
    my @q1;
    my @q2;
    while( <$fh> ) {
        chomp;
        if(/estimated ln prob of data/i) {
            $prob->[0] = [split("=", $_)]->[1];
        } elsif(/mean value of ln likelihood/i) {
            $prob->[1] = [split("=", $_)]->[1];
        } elsif(/variance of ln likelihood/i) {
            $prob->[2] = [split("=", $_)]->[1];
        } elsif(/mean value of alpha/i) {
            $prob->[3] = [split("=", $_)]->[1];
        } elsif(/proportion of membership of/i) {
            for my $i (1..5) {
                my $tmp = readline($fh);
            }
            while( my $line = readline($fh) ) {
                chomp($line);
                last if $line =~ /^[\-\s]+$/;
                push @q1, $line;
            }
        } elsif(/inferred ancestry of individuals/i) {
            for my $i (1..1) {
                my $tmp = readline($fh);
            }
            while( my $line = readline($fh) ) {
                chomp($line);
                last if $line =~ /^[\-\s]*$/;
                push @q2, $line;
            }
        }
    }
    return ($prob, \@q1, \@q2);
}
sub plot_structure {
    my ($p) = @_;
    my ($tag, $nl, $np, $ni) = map {$p->{$_}} qw/tag nl np ni/;
    my $dirE = dir($DIR_Src, "distruct1.1");
    chdir($dirE) or die "cannot cd to $dirE\n";
    my $dirW = dir($DIR_Misc3, "structure", $tag);
    my $d04 = file($dirW, "04_qmatrix");
    my $d05 = file($dirW, "05_stat");
    my $d06 = file($dirW, "06_figs");
    system("mkdir -p $d06") unless -d $d06;
    for my $k (@{$p->{ks}}) {
        for my $i (@{$p->{reps}}) { 
            my $f041 = file($d04, "k_$k\_rep_$i\_pop.txt");
            my $f042 = file($d04, "k_$k\_rep_$i\_ind.txt");
            next unless -s $f041;
            system("cp $f041 tmp_pop.txt");
            system("cp $f042 tmp_ind.txt");
            my $cmd = "distruct -K $k -M $np -N $ni -p tmp_pop.txt -i tmp_ind.txt -b none -a none -c none";
            open(JJ, $cmd." |") || die "Failed: $! in \n$cmd\n";
            while ( <JJ> ){
                chomp;
                print $_."\n";
            }
            my $f05 = file($d05, "k_$k\_rep_$i.ps");
            system("mv casia.ps $f05");
            print join("\t", "K=$k", "Rep=$i", "Okay")."\n";
        }
    }
}


