#!/bin/bash

cd $data/genome/HM056


awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR!=1) {$17="gene";} print}' 23.gtb > ../40_gene.gtb
cd ..

awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR!=1) {$18=$17; $17="CRP"; print}}' spada.crp.gtb | cut -f1-18 > 43_crp.gtb

gtb2bed.pl -i 42_nbs.gtb -o 42_nbs.bed
gtb2bed.pl -i 43_crp.gtb -o 43_crp.bed

cat 4[03]*.gtb > 49.gtb
gtbdedup.pl -i 49.gtb -o 51.gtb


gtb2gff.pl -i 51.gtb -o 51.gff
gtb2bb.pl -i 51.gtb -s 15.sizes -o 51.bb
gtb2bed.pl -i 51.gtb -o 51.bed
gtb2fas.pl -d 11_genome.fas -i 51.gtb -o 51.fas
#gtb2tbl.pl -i 51.gtb -o 51.tbl
