#!/usr/bin/perl -w
unshift (@INC, "INC");
require "my.pl";
use strict;
use Switch;
our $dbh;
#将group_id_arr转换为snp_id_arr
sub group_to_snp() {
    my ($ref) = @_;
    our $output_mode;
    my @snp_id_arr;
    foreach my $group_id (@$ref) {
        if($output_mode==1) {
            my $query = "SELECT b.id FROM mt_seq AS a, mt_snp AS b"
                ." WHERE a.snp_id=b.id AND a.group_id=$group_id";
            my $sth = $dbh->prepare($query);
            $sth->execute();
            while(my $rs = $sth->fetchrow_array()) {
                push (@snp_id_arr,$rs);
            }
        } else {
            my $query = "SELECT * FROM mt_group WHERE id=$group_id";
            my $sth = $dbh->prepare($query);
            $sth->execute();
            my $hashref = $sth->fetchrow_hashref();
            my @id_str_arr = split(",",$hashref->{"snp_id"});
            foreach my $id_str (@id_str_arr) {
                my @id_pair = split("-",$id_str);
                if(@id_pair == 1) {
                    push (@snp_id_arr,$id_pair[0]);
                } else {
                    for(my $j=$id_pair[0];$j<=$id_pair[1];$j++) {
                        push (@snp_id_arr,$j);
                    }
                }
            }
        }
    }
    return @snp_id_arr;
}
#对("1=4","1=7","2=3","3=6","4=6")返回("1=4=7","2=3=4=6")
sub cluster_id_pair() {
    my ($ref) = @_;
    my @group_arr;
    for(my $i=0; $i<@$ref; $i++) {
        my ($g,$h) = split("=",$ref->[$i]);
        my $flag = 0;
        for(my $j=0; $j<@group_arr; $j++) {
            if($group_arr[$j]=~/ $g / || $group_arr[$j]=~/ $h /) {
                if($group_arr[$j]!~/ $g /) {
                    $group_arr[$j] .= "$g ";
                }
                if($group_arr[$j]!~/ $h /) {
                    $group_arr[$j] .= "$h ";
                }
                $flag = 1;
            }
        }
        if($flag == 0) {
            push(@group_arr, " $g $h ");
        }
    }
    return @group_arr;
}
#对两个以空格键分格的字符串返回相差的个数
sub compare() {
    my ($string1, $string2) = @_;
    my @ary1 = split(" ",$string1);
    my @ary2 = split(" ",$string2);
    my $count=0;
    foreach my $variant (@ary1) {
        if($string2 !~ /$variant/i) {
            $count ++;
        }
    }
    foreach my $variant (@ary2) {
        if($string1 !~ /$variant/i) {
            $count ++;
        }
    }
    return $count;
}
sub get_snp_by_id() {
    my ($snp_id) = @_;
    my $query = "SELECT snp FROM mt_snp WHERE id=$snp_id";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $rs = $sth->fetchrow_array();
    return $rs;
}
#对输入的snp字符串进行range和ambiguous characters的check并返回重构的字符串
sub reconstruct_snp_string() {
    my ($snp_str) = @_;
    my @snp_arr = split(" ",$snp_str);
    my $buffer = "";
    foreach my $snp (@snp_arr) {
        if($snp =~ /^([.0-9]+)([A-Zd])$/) {
            $buffer .= $1.$2." " if &check_ele($1,$2)==1;
        } else {
            $buffer = "error snp string value";
        }
    }
    return $buffer; 
}
#对输入的单个snp字符串进行range和ambiguous characters的check
sub check_ele() {
    my ($pos, $value) = @_;
    our $range_include_string;
    our $range_exclude_string;
    our $ignore_chr_string;
    #check position
    my @range_included = split(",",$range_include_string);
    my @include_begin;
        my @include_end;
        if(@range_included) {
        foreach my $ele (@range_included) {
            my @pair = split("-",$ele);
            push (@include_begin,$pair[0]);
            push (@include_end,$pair[1]);
        }
    } else {
        push (@include_begin,1);
        push (@include_end,16569);
    }
    my @range_excluded = split(",",$range_exclude_string);
    my @exclude_begin;
    my @exclude_end;
    if(@range_excluded) {
        foreach my $ele (@range_excluded) {
            my @pair = split("-",$ele);
            push (@exclude_begin,$pair[0]);
            push (@exclude_end,$pair[1]);
        }
    } else {
        push (@exclude_begin,0);
        push (@exclude_end,0);
    }
    
    my $flag_pos = 0;
    my $i = 0;
    foreach (@include_begin) {
        if($pos>=$_ && $pos<=$include_end[$i]) {
            $flag_pos = 1;
        }
        $i++;
    }
    my $j = 0;
    foreach (@exclude_begin) {
        if(int($pos)>=$_ && int($pos)<=$exclude_end[$j]) {
            $flag_pos = 0;
        }
        $j++;
    }
    
    my $flag_value = 1;
    my @ignore_chrs = split(",",$ignore_chr_string);
    foreach my $ignore_chr (@ignore_chrs) {
        if($value =~ /$ignore_chr/) {
            $flag_value = 0;
        }
    }
    
    if($flag_pos==1 && $flag_value==1) {
        return 1;
    } else {
        return 0;
    }
}
sub base_transition() {
    my ($nt) = @_;
    switch ($nt) {
        case "A" { return "G"; }
        case "G" { return "A"; }
        case "T" { return "C"; }
        case "C" { return "T"; }
        else { return "?"; }
    }
}
#删除数组中重复元素
sub delete_repeat() {
    my ($ref) = @_;
    my @array_clear;
    my $i=0;
    my $tmp = " ";
    while($ref->[$i]) {
        if($tmp !~ / $ref->[$i] /i) {
            $tmp .= $ref->[$i]." ";
            push (@array_clear,$ref->[$i]);
        }
        $i ++;
    }
    return @array_clear;
}
#根据诸如"DQ11111-DQ11122,DQ22222,DQ22233-DQ33333"形式的Accession Number字符串重构一个array
sub acc_expand() {
    my ($acc) = @_;
    my @acc_str_arr = split(",",$acc);
    my @acc_arr;
    foreach my $acc_str (@acc_str_arr) {
        my @pair = split("-",$acc_str);
        if(@pair==1) {
            push (@acc_arr,$acc_str);
        } elsif(@pair==2) {
            $pair[0] =~ /\b([A-Z]{1,3})(\d+\b)/i;
            my $prefix = $1;
            my $start = $2;
            my $surfix_len = length($start);
            $pair[1] =~ /\b([A-Z]{1,3})(\d+\b)/i;
            if($1 ne $prefix) {
                die("Error in processing ".$acc_str."\n");
            }
            my $stop = $2;
            if($start>$stop) {
                my $tmp = $start;
                $start = $stop;
                $stop = $tmp;
            }
            for(my $i=$start;$i<=$stop;$i++) {
                push (@acc_arr, $prefix.sprintf("%0".$surfix_len."d",$i));
            }
        } else {
            die("Error in processing ".$acc_str."\n");
        }
    }
    return @acc_arr;
}
#将形如"1-4,6,8-9"的字符串转换为数组(1,2,3,4,6,8,9)
sub expand() {
    my ($range_str) = @_;
    my @range_pair_arr = split(",",$range_str);
    my @id_arr;
    foreach my $range_pair (@range_pair_arr) {
        my @pair = split("-",$range_pair);
        if(@pair==1) {
            push (@id_arr,$range_pair);
        } elsif(@pair==2) {
            my $start = $pair[0];
            my $stop = $pair[1];
            if($start>$stop) {
                my $tmp = $start;
                $start = $stop;
                $stop = $tmp;
            }
            for(my $i=$start;$i<=$stop;$i++) {
                push (@id_arr, $i);
            }
        } else {
            die("Error in processing ".$range_pair."\n");
        }
    }
    return @id_arr;
}
#将形如(1,2,3,4,6,8,9)的数组转换为字符串"1-4,6,8-9"
sub contract() {
    my ($ref) = @_;
    my @range_arr;
    my $begin = -3;
    my $end = -2;
    foreach my $i (sort {$a<=>$b} @$ref) {
        if($i == $end+1) {
            $end = $i;
        } else {
            if($begin>0) {
                if($begin==$end) {
                    push (@range_arr,$begin);
                } else {
                    push (@range_arr,$begin."-".$end);
                }
            }
            $begin = $i;
            $end = $i;
        }
    }
    if($begin>0) {
        if($begin==$end) {
            push (@range_arr,$begin);
        } else {
            push (@range_arr,$begin."-".$end);
        }
    }
    return join(",",@range_arr);
}
#为某个snp分配一个唯一的位置
my %snp_pos_arr;
sub get_pos() {
    my ($snp) = @_;
    if($snp =~ /^([.\d]+)([A-Zd])/) {
        my ($pos, $value) = ($1, $2);
        my $i=0;
        while(exists($snp_pos_arr{$pos+$i*0.01}) && $snp_pos_arr{$pos+$i*0.01} ne $snp) {
            $i++;
        }
        if(!exists($snp_pos_arr{$pos+$i*0.01})) {
            $snp_pos_arr{$pos+$i*0.01} = $snp;
        }
    } else {
        die("Illegal SNP string value\n");
    }
}

return 1;
