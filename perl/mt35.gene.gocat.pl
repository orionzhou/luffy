#!/usr/bin/perl -w
use strict; 
use Common; 
use GO::Parser;
use Parser;
use Path::Class;
use Data::Dumper;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

#updateGoFiles();
my $dir = "$ENV{'misc1'}/cat/mt_35";

my $f01 = "$dir/01_gene_ipr.tbl";
#get_ipr_ids($f01);
my $f02 = "$dir/02_gene_gos.tbl";
my $f03 = "$dir/03_gene_go.tbl";
#map_ipr_to_go($dir, $f01, $f02, $f03);

sub get_ipr_ids {
  my ($fo) = @_;
  my $fi = "$ENV{'genome'}/Mtruncatula_35/53.gtb";
  my $t = readTable(-in=>$fi, -header=>1);
  open(my $fh, ">$fo") or die "cannot write $fo\n";;
  print $fh join("\t", qw/gene_id interpro_ids/)."\n";
  for my $i (0..$t->nofRow-1) {
    my ($idM, $note) = map {$t->elm($i, $_)} qw/idM note/;
    my @iprs;
    while($note =~ /(IPR[\d]+)/g) {
      push @iprs, $1;
    }
    print $fh join("\t", $idM, join(" ", @iprs))."\n";
  }
}
sub map_ipr_to_go {
  my ($dir, $fi, $fo1, $fo2) = @_;
  my $f0 = "$dir/../go/interpro2go";
  my $h = parse_interpro2go($f0);
  
  my $t = readTable(-in=>$fi, -header=>1);
  my $fh1 = new IO::File $fo1, "w";
  print $fh1 join("\t", qw/gene_id go_ids/)."\n";
  my $fh2 = new IO::File $fo2, "w";
  print $fh2 join("\t", qw/gene_id go_id/)."\n";
  for my $i (0..$t->nofRow-1) {
      my ($id_gene, $ids_ipr) = $t->row($i);
      next unless $ids_ipr;
      my @ids_ipr = split(" ", $ids_ipr);
      my @ids_go;
      for my $id_ipr (@ids_ipr) {
          next unless exists $h->{$id_ipr};
          push @ids_go, @{$h->{$id_ipr}};
      }
      @ids_go = uniq(@ids_go);
      next unless @ids_go;
      print $fh1 join("\t", $id_gene, join(" ", @ids_go))."\n";
      print $fh2 join("\n", map {join("\t", $id_gene, $_)} @ids_go)."\n";
  }
  close $fh1;
  close $fh2;
}
  sub parse_interpro2go {
    my ($fi) = @_;
    my $fh = new IO::File $fi, "r";
    my $h;
    while(<$fh>) {
        chomp;
        next if /^\!/;
        if(/(IPR\d+).+(GO\:\d+)/) {
            $h->{$1} = [] unless exists $h->{$1};
            push @{$h->{$1}}, $2;
        } else {
            die "cannot parse:\n".$_."\n";
        }
    }
    close $fh;
    return $h;
} 

$f02 = file($dirW, "02_desc.txt");
#getFeDetail(-out=>$f02, -refdb=>$refDb, -ids=>\@ids);
my $f11 = file($dirW, "11_cat_ori.txt");
my $f21 = file($dirW, "21_cat_antoine.txt");
#catRefine($f11, $f21);
my $f22 = file($dirW, "22_cat_go.txt");
my $f23 = file($dirW, "23_cat_nbs.txt");
#catNbs($f23);
my $f24 = file($dirW, "24_cat_crp.txt");
#catCrp($f24);
my $f31 = file($dirW, "31_cat_info.txt");
my $f32 = file($dirW, "32_cat.txt");
#cat1(-ins=>[$f21, $f22, $f23, $f24], -out=>$f31);
#cat2(-ins=>[$f21, $f22, $f23, $f24], -ids=>\@ids, -cat=>$f31, -out=>$f32);
sub cat1 {
    my ($fis, $fo) = rearrange(['ins', 'out'], @_);
    my $catH;
    for my $i (0..@$fis-1) {
        my $fi = $fis->[$i];
        my $t = readTable(-in=>$fi, -header=>1);
        for my $catStr ($t->col("family")) {
            my @cats = split(";", $catStr);
            for my $cat (@cats) {
                $catH->{$cat} ||= [$i+1, 0];
                $catH->{$cat}->[1] ++;
            }
        }
    }
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/familyId family dataset nGene/)."\n";
    my $cnt = 0;
    for my $cat (sort keys %$catH) {
        print $fh join("\t", ++$cnt, $cat, @{$catH->{$cat}})."\n";
    }
}
sub cat2 {
    my ($ids, $fis, $fc, $fo) = rearrange(['ids', 'ins', 'cat', 'out'], @_);
    my $tc = readTable(-in=>$fc, -header=>1);
    my $catH = {map {$tc->elm($_, "family") => $tc->elm($_, "familyId")} (0..$tc->nofRow-1) };
    my $h;
    for my $fi (@$fis) {
        my $t = readTable(-in=>$fi, -header=>1);
        for my $i (0..$t->nofRow-1) {
            next unless $t->elm($i, "family");
            my $id = $t->elm($i, "id");
            $h->{$id} ||= [];
            my @cats = split(";", $t->elm($i, "family"));
            my @catIds;
            for my $cat (@cats) {
                die "$cat not in catH\n" unless exists $catH->{$cat};
                push @catIds, $catH->{$cat};
            }
            push @{$h->{$id}}, @catIds;
        }
    }
    my $fho = new IO::File $fo, "w";
    print $fho join("\t", qw/id family/)."\n";
    for my $id (@$ids) {
        my $cats = $h->{$id};
        if( $cats ) {
            for (@$cats) {
                print $fho join("\t", $id, $_)."\n";
            }
        }
    }
}
sub catRefine {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $fams = $t->delCol("family");
    $fams = [ map {s/(^\s+)|(\s+$)//g ? $_ : $_} @$fams ];
    $fams = [ map {s/\s*\;\s*/\;/g ? $_ : $_} @$fams ];
    $t->addCol($fams, "family");
    my @cats = uniq(map {split ";"} @$fams);
    for my $cat (sort @cats) {
        my @matches = grep /\Q$cat\E/i, @$fams;
        print join("\t", scalar(@matches), $cat)."\n";
    }
    $t->delCol("desc");
    $t->sort("family", 1, 0);
    my $fho = new IO::File $fo, "w";
    print $fho $t->tsv(1);
}
sub catNbs {
    my ($fo) = @_;
    my $fi = file($DIR_Misc2, "nbs", "mt_30", "01_cat.txt");
    my $t = readTable(-in=>$fi, -header=>1);
    my %h = map {$t->elm($_, "Gene model") => $t->elm($_, "cat")} (0..$t->nofRow-1);
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id family/)."\n";
    for my $str (keys %h) {
        my $cat = $h{$str};
        $str =~ s/\s+//g;
        my @ids = split(",", $str);
        for (@ids) {
            print $fh join("\t", $_, "NBS_$cat")."\n";
        }
    }
}
sub catCrp {
    my ($fo) = @_;
    my $fi = file($DIR_Misc2, "ncr", "mt_30", "02_hmm", "07_best_hits.txt");
    my $t = readTable(-in=>$fi, -header=>1);
    $t = $t->match_pattern("\$_->[6] < 1");
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/id family/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam) = map {$t->elm($i, $_)} qw/hId qId/;
        $fam =~ /(\d+)/;
        my $cat;
        if($1 <= 1030) {
            $cat = "CRP0000-1030";
        } elsif($1 >= 1040 && $1 <= 1530) {
            $cat = "CRP1040-1530";
        } else {
            die unless $1 <= 6250;
            $cat = "CRP1600-6250";
        }
        print $fh join("\t", $id, $cat)."\n";
    }
}
sub getAllIds {
    my $ld = Localdb->new(-db=>$refDb);
    my @fes = $ld->getFeatures(-types=>'mRNA');
    my $f01 = file($DIR_Misc1, "cat", "61_all.txt");
    my $fh = new IO::File $f01, "w";
    print $fh join("\t", qw/id chr start end evidence description/)."\n";
    for (sort {$a->seq_id cmp $b->seq_id || $a->start <=> $b->start} @fes) {
        next unless $_->seq_id =~ /^chr/i;
        print $fh join("\t", $_->id, $_->seq_id, $_->start, $_->end, $_->source, join(" ", $_->get_tag_values("Note")))."\n";
    }
}


sub updateGoFiles {
    my $dir = dir($DIR_Misc1, "cat", "go");
    my @urls = (
        "http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo", 
        "http://archive.geneontology.org/latest-termdb/go_daily-termdb.obo-xml.gz",
        "http://www.geneontology.org/external2go/interpro2go"
    );
    for my $url (@urls) {
        runCmd("wget --directory-prefix=$dir -N $url");
    }
    my $f = file($dir, "go_daily-termdb.obo-xml.gz");
    runCmd("gunzip -f $f");
}
sub go_classify {
    my ($fi, $fo1, $fo2) = @_;
    my $fGo = file($DIR_Misc1, "cat/go/go_daily-termdb.obo-xml");
    my $t = readTable(-in=>$fi, -header=>1);

    my $parser = GO::Parser->new({handler=>"obj", use_cache=>1});
    $parser->parse($fGo);
    my $graph = $parser->handler->graph;

    my ($ref1, $ref2, $catName) = ({}, {}, {});
    my $catLevel = 3;
    my $cnt = 0;
    for my $i (0..$t->nofRow-1) {
        my ($id, $go_ids) = map {$t->elm($i, $_)} qw/gene_id go_ids/;
        my @go_ids = split(" ", $go_ids);
        for my $go_id (@go_ids) {
            my $term = $graph->get_term($go_id);
            if(!$term) {
                print "$go_id not found\n";
                next;
            }
            $catName->{$term->acc} = $term->name if !exists($catName->{$term->acc});
#      print $fh "\t".join("\t", $term->acc."[".$term->name."]",$score)."\n";
#      my $parents = $term->get_parent_terms;
            my $paths = $graph->paths_to_top($go_id);
            for my $path (@$paths) {
                my $termAry = $path->term_list();
                my @tmpAry;
                for my $t (@$termAry) {
                    #push (@tmpAry, $t->acc."[".$t->name."]");
                    push (@tmpAry, $t->acc);
                    $catName->{$t->acc} = $t->name if !exists($catName->{$t->acc});
                }
                my @tmpAry2 = reverse @tmpAry;
                push (@tmpAry2, $term->acc);
                my $catLevelA = scalar(@tmpAry2)<$catLevel ? @tmpAry2 : $catLevel;
                my $cat = $tmpAry2[$catLevelA-1];
                $ref1->{$cat} = [] unless exists $ref1->{$cat};
                push @{$ref1->{$cat}}, $id;
                $ref2->{$id} = [] unless exists $ref2->{$id};
                push (@{$ref2->{$id}}, $cat);
            }
        }
        last if ++$cnt < 1;
    }
    my $fh1 = new IO::File $fo1, "w";
    print $fh1 join("\t", qw/gene_id go_ids/)."\n";
    for my $id (sort(keys(%$ref2))) {
        my @gos = uniq(@{$ref2->{$id}});
        print $fh1 join("\t", $id, join(" ", @gos))."\n";
    }
    my $fh2 = new IO::File $fo2, "w";
    print $fh2 join("\t", qw/gene_id go_id/)."\n";
    for my $go_id (sort(keys(%$ref1))) {
        my @ids = uniq(@{$ref1->{$go_id}});
        next unless @ids >= 45;
        for my $id (sort(@ids)) {
            print $fh2 join("\t", $id, $go_id, $catName->{$go_id})."\n";
        }
    }
}


my $f80 = file($DIR_Misc2, "version", "08_sum.txt");
my $f81 = file($dirW, "81_func.txt");
#testRun(-in1=>$f32, -in2=>$f80, -db=>"mt_35", -out=>$f81);
sub testRun {
    my ($fi1, $fi2, $db, $fo) = rearrange([qw/in1 in2 db out/], @_);
    my $ld = Localdb->new(-db=>$db);
    my $t1 = readTable(-in=>$fi1, -header=>1);
    my $t2 = readTable(-in=>$fi2, -header=>1);
    $t1 = $t1->match_pattern("\$_->[1] == 34");
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/mt_30_id mt_35_id mt_35_desc/)."\n";
    for my $i (0..$t1->nofRow-1) {
        my $id1 = $t1->elm($i, "id");
        $id1 =~ s/\.\d$//;
        my $t3 = $t2->match_pattern("\$_->[0] eq '$id1'");
        die "$id1 not found\n" unless $t3->nofRow;
        for my $j (0..$t3->nofRow-1) {
            unless($t3->elm($j, "mt_35_id")) {
                print $fh join("\t", $id1, "", "")."\n";
                next;
            }
            my $idStr = $t3->elm($j, "mt_35_id");
            my @descs;
            for my $str (split(" ", $idStr)) {
                $str =~ /^(\w+)/;
                my $fe = $ld->getFeatureByName("$1.1");
                push @descs, join(" ", $fe->get_tag_values("Note"));
            }
            print $fh join("\t", $id1, $idStr, join("| ", @descs))."\n";
        }
    }
}


sub filter_mt_30_goterms {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my @ids = $t->col("id");
    print @ids." items\n";
    @ids = uniq(@ids);
    print @ids." gene ids\n";
    @ids = grep /^Medtr.*/, @ids;
    print @ids." chr1-8 gene ids\n";
    my ($cnt1, $cnt2) = (0) x 2;
    my $scoreCutOff = 100;
    my $fh = new IO::File $fo, "w";
    print $fh join("\t", qw/gene_id go_ids/)."\n";
    for my $id (@ids) {
        my $t2 = $t->match_pattern("\$_->[0] eq '$id' && \$_->[3] > $scoreCutOff");
        if( $t2->nofRow > 0 ) {
            $cnt1 ++;
            $cnt2 += $t2->nofRow;
            my @groups = $t2->col("group");
            print $fh join("\t", $id, join(" ", @groups))."\n";
        }
    }
    print "$cnt1 genes with $cnt2 go items with score >= $scoreCutOff\n";
}



