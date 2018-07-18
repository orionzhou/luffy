#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
use Common; 
use Seq;
use Medicago;
use Data::Dumper;

my $dir = "/home/youngn/zhoup/Data/misc1/seq06";
my $f01 = "$dir/01_id.tbl";

my $f00 = "/home/youngn/zhoup/Data/genome/Mtruncatula_3.5/21_gene.gtb";
#awk 'BEGIN{FS="\t"; print "id"} {if($3 ~ /chr[1-8]/ && $11 !~ /transposable_element_gene/) {print $2}}' $f00 > $f01

#seqRec -i 01_id.tbl -o 13_seq -t acc84 -c 1
#print "awk '{ if( NR==1 || (NR>1 && $1!=\"id\") ) print}' $dir/13_seq_chr[1-8].tbl > $dir/13_seq.tbl\n";
#print "rm $dir/13_seq_chr[1-8].tbl\n\n";
#print "cat $dir/13_seq_chr[1-8].log > $dir/13_seq.log\n";
#print "rm $dir/13_seq_chr[1-8].log\n\n";

#cut -f1-6 13_seq.tbl > 32_info.tbl

my $opt_ind = "acc31";
my ($accs, $opt_conf) = get_mt_ids($opt_ind);
print "seqSplit -i 13_seq.tbl -t $opt_ind -c $opt_conf -d $opt_ind/01_stat.tbl -o $opt_ind\n"; 

my ($sub_pre, $sub_opt1, $sub_opt2, $sub_opt3) = ("subset_opt12", "opt12.1", "opt12.2", "opt12.3");
#subset_fasta_by_opt(-dir=>"$dir/$opt_ind", -pre=>$sub_pre, -opts=>[$sub_opt1, $sub_opt2, $sub_opt3]);
sub subset_fasta_by_opt {
    my ($dir, $pre, $opts) = rearrange(['dir', 'pre', 'opts'], @_);
    system("mkdir -p $dir/$pre") unless -d "$dir/$pre";
    my $f_id = "$dir/02_ids.tbl";
    my $ids = getIds($f_id);

    my @types = qw/cds intron utr3 utr5/;
    my $h_opt;
    for my $opt (@$opts) {
        for my $type (@types) {
            my $subdir = "$dir/$pre/$opt/$type";
            system("mkdir -p $subdir") unless -d $subdir;
        }
        $h_opt->{$opt} = [get_mt_ids($opt)]->[0];
    }

    for my $i (0..@$ids-1) {
        my $id = $ids->[$i];
        for my $type (@types) {
            my $fi = "$dir/$type/$id.fas";
            next unless -s $fi;
            my $h_seq = readSeq($fi, 2);
            for my $opt (@$opts) {
                my $fo = "$dir/$pre/$opt/$type/$id.fas";
                my $accs = $h_opt->{$opt};
                my @seqs = map { Bio::Seq->new(-id=>$_, -seq=>$h_seq->{$_}) } @$accs;
                writeSeq(\@seqs, $fo);
            }
        }
        printf " %5d / %5d done\r", $i+1, scalar(@$ids);
        last if $i < 0;
    }
    print "\n";
}

