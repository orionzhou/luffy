
ln -sf augustus/41.gtb 41.gtb
awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR>1) {$15="mRNA"; $16="nbs"; print}}' 42.nbs/11.gtb | cut -f1-18 > 42.nbs.gtb
awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR>1) {$15="mRNA"; $16="crp"; $18=""; print}}' spada.crp.gtb | cut -f1-18 > 43.crp.gtb

cat 4[1-3]*.gtb > 49.gtb
gtb.dedup.pl -i 49.gtb -o 51.gtb

awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR==1 || tolower($16) != "te") print}' 51.gtb > 55_noTE.gtb

gtb2gff.pl -i 51.gtb -o 51.gff
gtb.idx.pl -i 51.gtb -s 15.sizes