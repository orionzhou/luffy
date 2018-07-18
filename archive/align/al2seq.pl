#!/usr/bin/perl -w
#use the bl2seq executable, align the rCRS with input sequence(Bio::Seq Object)
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Tools::Run::StandAloneBlast;

my $working_dir = "./";
my $seq_path = "in";
my $result_path = "out";
my $rCRS_path = "rCRS.gb";

my $rCRS_file = Bio::SeqIO->new(-file => $working_dir.$rCRS_path, -format=>'genbank');
my $rCRS = $rCRS_file->next_seq();
my $seq_file = Bio::SeqIO->new(-file=>$working_dir.$seq_path, -format=>'fasta')
	or die("Cannot open sequence file");
open(OUT, ">".$result_path) or die("cannot write to outfile");
while(my $seq = $seq_file->next_seq()) {
	&al2seq($seq, $rCRS);
}

sub al2seq() {
    my ($seq, $rCRS) = @_;
    my $align_result_file = $working_dir."bl2seq.out";
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
              $mut_arr{$loc_in_rCRS->start()} = uc(substr($seq->seq(),$loc_in_seq->start()-1,1));
            } elsif($loc_in_seq->length() == 0) {
              #print "deletion ".$loc_in_rCRS->start()." -> [".$loc_in_seq->start()."_".$loc_in_seq->end()."]";
              $mut_arr{$loc_in_rCRS->start()} = "d";
            } else {
              print $loc_in_seq->length()." column not in query seq";
              exit;
            }
          } elsif($loc_in_rCRS->length() == 0) {
            #print "insertion [".$loc_in_rCRS->start()."_".$loc_in_rCRS->end()."] -> ".$loc_in_seq->start();
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
        print OUT $rCRS_aligned->start()."\t".$rCRS_aligned->end()."\t";
      }
    }
    foreach (sort {$a <=> $b } keys(%mut_arr)) {
		print OUT $_."$mut_arr{$_}\t";
    }
    print OUT "\n";
}