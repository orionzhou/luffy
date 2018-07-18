#!/bin/bash

cd $data/genome/Mtruncatula_4.0/raw

ln -sf JCVI.Medtr.v4.20130313.fasta 01_all.fas
#grep '>' 01_all.fas | grep -o 'chr[0-9]\+' > 02.chrs.txt
#grep '>' 01_all.fas | grep -o 'scaffold[0-9]\+' > 02.scaffolds.txt
#seqret.pl -d 01_all.fas -b 02.chrs.txt -o 03.chrs.fas
#seqret.pl -d 01_all.fas -b 02.scaffolds.txt -o 03.scaffolds.fas
#seqconcat.pl -i 03.scaffolds.fas -o 03.scaffolds.concat.fas -p 03.scaffolds.concat.tbl -d chrU -g 1000

#cat 03.chrs.fas 03.scaffolds.concat.fas > ../11_genome.fas
seq.mt40.pl -i 01_all.fas -o 05.fas -p 05.gal
gal2psl.pl -i 05.gal -o 05.psl
pslToChain 05.psl 05.chain

gff.jcvi.pl -i 11.gff -o 12_jcvi_fixed.gff
gff.coord.pl -i 12_jcvi_fixed.gff -p 03.scaffolds.concat.tbl -o 15_global_loc.gff
gff2gtb.pl -i 15_global_loc.gff -o 15_global_loc.gtb
gtb.phase.pl -i 15_global_loc.gtb -s ../11_genome.fas -o 17_phase_fixed.gtb
gtb.dedup.pl -i 17_phase_fixed.gtb -o 25.dedup.gtb
#gtb.pickalt.pl -i 25.dedup.gtb -o 26.longest.gtb

ln -sf 31.gtb 25.dedup.gtb
gtb2gff.pl -i 31.gtb -o 31.gff
gtb2fas.pl -i 31.gtb -o 31.fas -d ../11_genome.fas


#run pipe_pfam() in mt.augus.pl -g HM101 