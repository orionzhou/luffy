require(plyr)
require(ape)
require(ggplot2)

options(stringsAsFactors = FALSE)

dirw = file.path(Sys.getenv("misc2"), "devel.atlas")

dirg = file.path(Sys.getenv("genome"), "Zmays_v4")
fg = file.path(dirg, "51.gtb")
#tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,16:18)]
#gb = group_by(tg, par)
#tg2 = summarise(gb, fam = names(sort(table(cat3), decreasing = T))[1])

### get data in
fi = file.path(dirw, "per_tissue_ave_RPM.txt")
ti = read.table(fi, sep = " ", header = T, quote = '"')

fn = file.path(dirw, "gnames.txt")
gnames = read.table(fn, header = T)[,1]

stopifnot(nrow(ti) == length(gnames))
tmap = data.frame(tname = rownames(ti), gname = gnames, stringsAsFactors = F)
fm = sprintf("%s/00.map.tsv", dirw)
write.table(tmap, fm, col.names = T, row.names = F, sep = "\t", quote = F)

rownames(ti) = gnames
fo = sprintf("%s/01.tsv", dirw)
write.table(ti, fo, col.names = T, row.names = T, sep = "\t", quote = F)

### run camoco on a lab workstation
source activate camoco
camoco --help
camoco build-refgen $genome/Zmays_v4/51.gff maize v34 v34 maize
camoco build-cob --rawtype RNASEQ --index-col 1 --min-single-sample-expr 0 --max-val 500 $misc2/devel.atlas/01.tsv sarah.da 1.0 maize
camoco health --out $misc2/devel.atlas/51.camoco/da sarah.da

source activate camoco
ipython
import camoco as co
x = co.COB("sarah.da")

x.clusters.to_csv("/home/springer/zhoux379/data/misc2/devel.atlas/51.clusters.csv")
#x.coex.score.to_csv(file.path(dirw, "51.camoco/coex.csv", index = False, float_format='%g')

### post-process cluster assignments
fc = file.path(dirw, "51.clusters.csv")
tc = read.csv(fc)

fm = sprintf("%s/00.map.tsv", dirw)
tm = read.table(fm, sep = "\t", header = T, stringsAsFactors = F)
nmap = tm$tname
names(nmap) = toupper(tm$gname)

tc = cbind(tc, tid = nmap[tc$X])
tc$cluster = tc$cluster + 1
fo = sprintf("%s/51.camoco/clusters.tsv", dirw)
write.table(tc[,-1], fo, col.names = T, row.names = F, sep = "\t", quote = F)
