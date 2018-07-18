#!/usr/bin/perl -w
#use the bl2seq executable, align the rCRS with input sequence(Bio::Seq Object)
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Tools::Run::StandAloneBlast;

my $kkk="jjj";
sub al2seq() {
    my ($seq, $rCRS) = @_;
    my $buffer = "";
    my $align_result_file = "./bl2seq.out";
    my $factory = Bio::Tools::Run::StandAloneBlast->new( -program => 'blastn',
    	-outfile => $align_result_file, -F => "F");
    my $bl2seq_report = $factory->bl2seq($seq, $rCRS);

    # Use AlignIO.pm to create a SimpleAlign object from the bl2seq report
    my $str = Bio::AlignIO->new(-file => $align_result_file, -format => 'bl2seq');
    my %mut_arr;
    while( my $aln = $str->next_aln() ) {
      if($aln->length > 100) {
        my $seq_aligned = $aln->get_seq_by_pos(1);
        my $rCRS_aligned = $aln->get_seq_by_pos(2);
        my $match_line = $aln->match_line();
        while($match_line =~ /([^\*])/g ) {
          my $unmatch_column = pos($match_line);
          my $loc_in_rCRS = $rCRS_aligned->location_from_column($unmatch_column);
          my $loc_in_seq = $seq_aligned->location_from_column($unmatch_column);
          if($loc_in_rCRS->length() == 1) {
            if($loc_in_seq->length() == 1) {
              $mut_arr{$loc_in_rCRS->start()} = uc(substr($seq->seq(),$loc_in_seq->start()-1,1));
            } elsif($loc_in_seq->length() == 0) {
              $mut_arr{$loc_in_rCRS->start()} = "d";
            } else {
              print $loc_in_seq->length()." column not in query seq";
              exit;
            }
          } elsif($loc_in_rCRS->length() == 0) {
            my @keys_mut = keys(%mut_arr);
            my $ins_num = 1;
            my $tmpstr = $loc_in_rCRS->start();
            foreach my $key_mut (@keys_mut) {
              if($key_mut =~ /^$tmpstr\./) {
                $ins_num ++;
              }
            }
            $mut_arr{$loc_in_rCRS->start().".".$ins_num} = substr($seq->seq(),$loc_in_seq->start()-1,1);
          } else {
            print "column not in rCRS seq";
            exit;
          }
        }
        $buffer .= $rCRS_aligned->start()."\t".$rCRS_aligned->end()."\t";
      }
    }
    foreach (sort {$a<=>$b} keys(%mut_arr)) {
		$buffer .= $_."$mut_arr{$_} ";
    }
    return $buffer;
}