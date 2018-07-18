#!/bin/bash

cd $data/genome/Mtruncatula_3.5/sequence

changeSeqID.pl -i 01_chr1-8.fa -o 01_chr1-8.renamed.fa

changeSeqID.pl -i 02_chr0.fa -o 02_chr0.renamed.fa
seqconcat.pl -i 02_chr0.renamed.fa -o 02_chr0.renamed.concat.fa -p 02_chr0.renamed.concat.tbl -gap 50000 -id chr0

changeSeqID.pl -i 03_chrU.fa -o 03_chrU.renamed.fa
seqconcat.pl -i 03_chrU.renamed.fa -o 03_chrU.renamed.concat.fa -p 03_chrU.renamed.concat.tbl -gap 1000 -id chrU

changeSeqID.pl -i 04_chrT.fa -o 04_chrT.renamed.fa
seqconcat.pl -i 04_chrT.renamed.fa -o 04_chrT.renamed.concat.fa -p 04_chrT.renamed.concat.tbl -gap 1000 -id chrT

cat 01_chr1-8.renamed.fa 02_chr0.renamed.concat.fa 03_chrU.renamed.concat.fa 04_chrT.renamed.concat.fa > ../11_genome.fa
cattab.pl -in 02_chr0.renamed.concat.tbl -in 03_chrU.renamed.concat.tbl -in 04_chrT.renamed.concat.tbl -out 12_assembly.tbl

cd ../gene

jcvi.pl -i 11.gff -o 12_jcvi.fixed.gff
gff_convert_loc.pl -i 12_jcvi_fixed.gff -p ../sequence/12_assembly.tbl -o 15_global_loc.gff
gff_to_gtb.pl -i 15_global_loc.gff -o ../21_gene.gtb -s ../11_genome.fa

cd ..
gtb_conv.pl -i 21_gene.gtb -o 21_gene.gff
gtb_conv.pl -i 21_gene.gtb -o 21_gene.fas -f seq2 -s ../11_genome.fa
gtb_conv.pl -i 21_gene.gtb -o 21_gene.tbl -f tbl