#!/usr/bin/perl -w
#use the bl2seq executable, align the rCRS with input sequence(Bio::Seq Object)
use strict;
use Bio::Seq;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Tools::Run::StandAloneBlast;

my $align_dir = "./";

sub al2seq() {
    my ($seq, $rCRS, @range) = @_;
    my @begin;
    my @end;
    if(@range) {
      foreach my $ele (@range) {
        my @pair = split("-",$ele);
        push (@begin,$pair[0]);
        push (@end,$pair[1]);
      }
    } else {
      push (@begin,1);
      push (@end,$rCRS->length);
    }
    #my $align_result_file = $align_dir.substr($oneseq,0,length($oneseq)-3).'.out';
    my $align_result_file = $align_dir."bl2seq.out";
    #print $rCRS->primary_id()."\t".$seq->primary_id()."\n";

    #call bl2seq executable
    my $factory = Bio::Tools::Run::StandAloneBlast->new( -program => 'blastn',
    	-outfile => $align_result_file, -F => "F");
    my $bl2seq_report = $factory->bl2seq($seq, $rCRS);

    # Use AlignIO.pm to create a SimpleAlign object from the bl2seq report
    my $str = Bio::AlignIO->new(-file => $align_result_file, -format => 'bl2seq');
    my %mut_arr;
    print $seq->id."\t";
    while( my $aln = $str->next_aln() ) {
      if($aln->length > 100) {
        #print $aln->length."\n";
        #print $aln->percentage_identity."\n";
        #print $aln->consensus_string()."\n";
        my $seq_aligned = $aln->get_seq_by_pos(1);
        my $rCRS_aligned = $aln->get_seq_by_pos(2);
#        my @params = (-refseq     => $rCRS->primary_id(),
#                           -allele1    => 'allele1',
#                           -allele2    => 'allele2',
#                           -delimiters => '[]',
#                           -separator  => '/');
#        print $aln->bracket_string(@params);

#        $aln->match_char(".");
#        $aln->missing_char("x");
#        $aln->gap_char("_");
        #print length($aln->match_line())."\n";
        my $match_line = $aln->match_line();
        #my $gap_line = $aln->gap_line();
        #print $match_line."\n".$gap_line."\n";
        while($match_line =~ /([^\*])/g ) {
          my $unmatch_column = pos($match_line);
          #print $unmatch_column.":";
          my $loc_in_rCRS = $rCRS_aligned->location_from_column($unmatch_column);
          my $loc_in_seq = $seq_aligned->location_from_column($unmatch_column);
          if($loc_in_rCRS->length() == 1) {
            if($loc_in_seq->length() == 1) {
              #print "mismatch ".$loc_in_rCRS->start()." -> ".$loc_in_seq->start();
              $mut_arr{$loc_in_rCRS->start()} = "[".substr($rCRS->seq(),$loc_in_rCRS->start()-1,1)
              	."/".substr($seq->seq(),$loc_in_seq->start()-1,1)."]";
            } elsif($loc_in_seq->length() == 0) {
              #print "deletion ".$loc_in_rCRS->start()." -> [".$loc_in_seq->start()."_".$loc_in_seq->end()."]";
              $mut_arr{$loc_in_rCRS->start()} = "[".substr($rCRS->seq(),$loc_in_rCRS->start()-1,1)
              	."/"."-]";
            } else {
              print $loc_in_seq->length()."column not in query seq";
              exit;
            }
          } elsif($loc_in_rCRS->length() == 0) {
            #print "insertion [".$loc_in_rCRS->start()."_".$loc_in_rCRS->end()."] -> ".$loc_in_seq->start();
            if( !exists $mut_arr{$loc_in_rCRS->start()."^".$loc_in_rCRS->end()} ) {
              $mut_arr{$loc_in_rCRS->start()."^".$loc_in_rCRS->end()} =
              	"[-/".substr($seq->seq(),$loc_in_seq->start()-1,1)."]";
            } else {
              my $tmp = substr($seq->seq(),$loc_in_seq->start()-1,1);
              $mut_arr{$loc_in_rCRS->start()."^".$loc_in_rCRS->end()} =~
              	s/\]/$tmp\]/;
            }
          } else {
            print "column not in rCRS seq";
            exit;
          }
        }
        print $rCRS_aligned->start()."\t".$rCRS_aligned->end()."\t";
      }
    }
        foreach (sort {$a <=> $b } keys(%mut_arr)) {
          $_ =~ /^\b(\d+)\b/;
          my $pos = $1;
          my $flag = 0;
          my $i = 0;
          foreach (@begin) {
            if($pos>=$_ && $pos<=$end[$i]) {
              $flag = 1;
            }
            $i++;
          }
          if($flag == 1) {
            print $_."$mut_arr{$_}\t";
          }
        }
    print "\n";
}
