awk 'BEGIN {FS="\t"; OFS="\t"} {print $6, $7-1, $8}' 23.gax | \
  sortBed -i stdin | mergeBed -i stdin > tmp23.bed
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2-1, $3}' 33.gax | \
  sortBed -i stdin | mergeBed -i stdin > tmp33.bed

awk 'BEGIN {FS="\t"; OFS="\t"} {print $6, $7-1, $8}' 31.3.gax | \
  sortBed -i stdin | mergeBed -i stdin > 31.3.bed
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2-1, $3}' 41.3.gax | \
  sortBed -i stdin | mergeBed -i stdin > 41.3.bed

subtractBed -a $genome/HM056/15.bed -b $genome/HM056/16_gap.bed | \
  subtractBed -a stdin -b tmp23.bed | \
  awk '($3-$2) >= 50' - > tmp.nov23.bed

subtractBed -a $genome/HM056/15.bed -b $genome/HM056/16_gap.bed | \
  subtractBed -a stdin -b 41.3.bed | bedfilter.pl -l 50 -o nov.bed


awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR != 1) print $1, $2-1, $3}' \
  nov3.tbl > nov3.bed