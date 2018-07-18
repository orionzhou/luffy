package SamPileup;
use strict; use Init; use Common; use Localdb; use Mapping; 
use Path::Class; use DBI; use Data::Dumper; use Statistics::Basic qw/:all/;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/bam2Vtb vtb2Fasta/;
@EXPORT_OK = qw//;
sub parseRegionString {
    my ($region) = @_;
    my ($chr, $beg, $end);
    if($region =~ /^([\d\w]+)\:([\d,]+)-([\d,]+)$/) {
        ($chr, $beg, $end) = ($1, $2, $3);
    } else {
        die "unknown region string: $region\n";
    }
    $beg =~ s/,//g;
    $end =~ s/,//g;
    return ($chr, $beg, $end);
}
sub bam2Vtb {
    my ($id, $region, $opt) = rearrange(['id', 'region', 'opt'], @_);
    my $dirW = dir($DIR_Repo, "mt_35");
    my $f_genome = file($dirW, "01_reference/41_genome.fa");
    my $f_bam = file($dirW, "10_bam", "05_pos_sorted", "$id.bam");

    my $cmd1 = "samtools mpileup -gf $f_genome $f_bam -r $region";
    my $cmd2 = "samtools mpileup -q1 -gf $f_genome $f_bam -r $region";
    my $cmd3 = "bcftools view -g -";
    my $lines1 = runCmd2("$cmd1 | $cmd3");
    my $lines2 = runCmd2("$cmd2 | $cmd3");
    my $f_tmp_1 = file($DIR_Tmp, "bampileup.vcf");
    write_file(-out=>$f_tmp_1, -ps=>[$lines1, $lines2], -opt=>1);

    my $t1 = vcf2Vtb($lines1, $region);
    my $t2 = vcf2Vtb($lines2, $region);
    my $f_tmp_2 = file($DIR_Tmp, "bamtable.tbl");
    write_file(-out=>$f_tmp_2, -ps=>[$t1, $t2], -opt=>2);

    return merge_vtb($t1, $t2);
}
sub write_file {
    my ($fo, $ps, $opt) = rearrange(['out', 'ps', 'opt'], @_);
    my $fh = new IO::File $fo, 'w';
    for my $p (@$ps) {
        if($opt == 1) {
            print $fh join("\n", @$p);
            print $fh "\n\n";
        } elsif($opt == 2) {
            print $fh $p->tsv(1);
            print $fh "\n";
        }
    }
    close $fh;
}
sub vcf2Vtb {
    my ($lines, $region) = @_;
    my ($chr, $beg, $end) = parseRegionString($region);

    my $co_mq = 10;
    my $co_cov_min = 2;
    my $co_cov_max = 100000;
    my $t = Data::Table->new([], [qw/chr pos ref alt indel type called cov quality/]);
    my %het = (AC=>'M', AG=>'R', AT=>'W', CA=>'M', CG=>'S', CT=>'Y', GA=>'R', GC=>'S', GT=>'K', TA=>'W', TC=>'Y', TG=>'K');

    my $last_pos = $beg;
    for (@$lines) {
        next if /^#/;
        my @t = split;
        die "not $chr: $t[0]\n" if ($chr ne $t[0]);
        die "[vcf2fq] unsorted input\n" if ($t[1] - $last_pos < 0);
        if ($t[1] - $last_pos > 1) {
            for my $pos ($last_pos..$t[1]-1) {
                $t->addRow( [$chr, $pos, '', '', '', 'nc', 'N', '', ''] );
            }
        }
        my ($pos, $ref, $alt) = @t[1,3,4];
        $alt = $1 if $alt =~ /^([A-Za-z.]+)(,[A-Za-z]+)*$/;
        my ($indel, $type, $called, $cov, $q, $mq, $fq) = ("") x 7; 
        $cov = $1 if $t[7] =~ /DP=(\d+)/;
        $mq = $1 if $t[7] =~ /MQ=(\d+)/;
        $fq = $1 if $t[7] =~ /FQ=(-?[\d\.]+)/;
        my $af = $1 if $t[7] =~ /AF1=([\d\.]+)/;
        $q = $fq;
        if (length($ref) == 1 && $t[7] !~ /INDEL/ && $t[4] =~ /^([A-Za-z.])(,[A-Za-z])*$/) {
            if ($q < 0) {
                my $tag = ($af < .5 || $alt eq ".");
                $called = $tag ? $ref : $alt;
                $type = "snp" if !$tag;
                $q = -$q;
            } else {
                $called = $het{"$ref$alt"};
                $called ||= 'N';
                $type = "het";
            }
        } elsif ($t[4] ne '.') {
            $indel = 1;
            if ($q < 0) {
                my $tag = ($af < .5 || $alt eq ".");
                $called = $tag ? $ref : $alt;
                if(!$tag) {
                    if(length($ref) > 1 && length($alt) == 1) {
                        $type = "del";
                    } elsif(length($ref) == 1 && length($alt) > 1) {
                        $type = "ins";
                    }
                }
                $q = -$q;
            } else {
                $called = 'N';
                $type = "het";
            }
        }
        $q = int($q + 33 + .499);
        $q = chr($q <= 126? $q : 126);
        unless ($cov >= $co_cov_min && $cov <= $co_cov_max && $mq >= $co_mq) {
            ($type, $called) = ("np", "N");
        }
        $t->addRow([$chr, $pos, $ref, $alt, $indel, $type, $called, $cov, $q]);
        $last_pos = $t[1];
    }
    if($last_pos < $end) {
        for my $pos ($last_pos+1..$end) {
            $t->addRow( [$chr, $pos, '', '', '', 'nc', 'N', '', ''] );
        }
    }
    return $t;
}
sub merge_vtb {
    my ($t1, $t2) = @_;
    $t2->header([qw/chr pos ref alt indel type called cov2 quality/]);
    $t2->addCol([('') x $t2->nofRow], 'cov1', 7);

    my $t1_info = { map { join(".", $t1->elm($_, "chr"), $t1->elm($_, "pos"), $t1->elm($_, "ref")) => $t1->rowRef($_) } (0..$t1->nofRow-1) }; 
    for my $i (0..$t2->nofRow-1) {
        my $key = join(".", map {$t2->elm($i, $_)} qw/chr pos ref/);
        my $row1 = $t1_info->{$key};
        my $cov1 = $row1->[7];
        $t2->setElm($i, "cov1", $cov1);
    }
    return $t2;
}
sub vtb2Fasta {
#chr pos ref alt indel type called cov1 cov2 quality
    my ($t) = @_;
    my $t1 = $t->match_pattern("\$_->[4] eq ''");
    my $t11 = $t1->match_pattern("\$_->[5] eq 'nc' || \$_->[5] eq 'np' || \$_->[5] eq 'het'");
    my $t12 = $t1->match_pattern("\$_->[5] eq 'snp'");
    my @snps = map {join(":", $t12->elm($_, "pos"), $t12->elm($_, "ref"), $t12->elm($_, "alt"))} (0..$t12->nofRow-1);
    my $t13 = $t1->match_pattern("\$_->[5] eq ''");
    my ($cnt_n, $cnt_s, $cnt_r) = ($t11->nofRow, $t12->nofRow, $t13->nofRow);
    die "$cnt_n + $cnt_s + $cnt_r != ".$t1->nofRow."\n" unless $cnt_n+$cnt_s+$cnt_r == $t1->nofRow;
    
    my $t2 = $t->match_pattern("\$_->[4] eq 1 && \$_->[5] eq 'ins'");
    my $cnt_i = $t2->nofRow;
    my @inss = map {join(":", $t2->elm($_, "pos"), $t2->elm($_, "ref"), $t2->elm($_, "alt"))} (0..$t2->nofRow-1);
    my $t3 = $t->match_pattern("\$_->[4] eq 1 && \$_->[5] eq 'del'");
    my $cnt_d = $t3->nofRow;
    my @dels = map {join(":", $t3->elm($_, "pos"), $t3->elm($_, "ref"), $t3->elm($_, "alt"))} (0..$t3->nofRow-1);

    my $beg = $t->elm(0, "pos");
    my @seqs;
    for my $i (0..$t1->nofRow-1) {
        my ($type, $nt) = map {$t1->elm($i, $_)} qw/type called/;
        $nt = "N" if $type eq "het";
        push @seqs, $nt;
    }
    my $seqStr1 = join("", @seqs);
    my $seqStr2 = applyIndel($beg, \@seqs, $t2, $t3);
    
    my $covs1 = [map {$_ ? $_ : 0} $t1->col("cov1")];
    my $covs2 = [map {$_ ? $_ : 0} $t1->col("cov2")];
    my $h = {
        n=>$cnt_n, snp=>$cnt_s, ref=>$cnt_r, ins=>$cnt_i, del=>$cnt_d,
        snps=>\@snps, inss=>\@inss, dels=>\@dels, 
        cov1_mean=>mean($covs1)->query, cov1_median=>median($covs1)->query, cov1_sd=>stddev($covs1)->query,
        cov2_mean=>mean($covs2)->query, cov2_median=>median($covs2)->query, cov2_sd=>stddev($covs2)->query,
    };

    return ($h, $seqStr1, $seqStr2);
}
sub applyIndel {
#pay attention to call insertion
    my ($beg, $seqs, $t1, $t2) = @_;
    for my $i (0..$t1->nofRow-1) {
        my ($pos, $ref, $alt, $type) = map {$t1->elm($i, $_)} qw/pos ref alt type/;
        my $idx = $pos - $beg;
        $seqs->[$idx] = $alt;
    }
    for my $i (0..$t2->nofRow-1) {
        my ($pos, $ref, $alt, $type) = map {$t2->elm($i, $_)} qw/pos ref alt type/;
        my $idx = $pos - $beg;
        for my $i ($idx+1..$idx+length($ref)-1) {
            $seqs->[$i] = "-";
        }
    }
    return join("", @$seqs);
}



1;
__END__
