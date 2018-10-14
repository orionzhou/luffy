require(plyr)
require(dplyr)
require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(xlsx)
require(RColorBrewer)
source("Location.R")
source("comp.fun.R")

dirw = file.path(Sys.getenv('misc3'), 'alpaca')
qnames_ac = c("HM034.AC", "HM056.AC", "HM340.AC")

##### work on individial Alpaca assemblies
org = "HM010"
orgp = "HM010"

dirg = file.path(Sys.getenv("genome"), org)

## check ovlp with segmental duplication blocks
fm = file.path(dirg, "raw.fix.fas.map")
tm = read.table(fm, header = F, sep = "\t", as.is = T)
colnames(tm) = c("ochr", "len", "nchr")

fg = file.path(dirg, "51.tbl")
tg = read.table(fg, header = F, sep = "\t", as.is = T)
colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
tgs = tg[tg$type == 'mrna',]
grg = with(tgs, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
grgn = grg[tgs$fam %in% c("TIR-NBS-LRR", "CC-NBS-LRR")]
grgc = grg[tgs$fam %in% c("CRP0000-1030", "NCR", "CRP1600-6250")]

f1 = sprintf("%s/nucmer_%s_alpaca_to_self/same_scaffold_repeat.1.bed", dirw, tolower(orgp))
t1 = read.table(f1, header = F, sep = " ", as.is = T)

f2 = sprintf("%s/nucmer_%s_alpaca_to_self/same_scaffold_repeat.2.bed", dirw, tolower(orgp))
t2 = read.table(f2, header = F, sep = " ", as.is = T)

sum(t1$V1 != t2$V1)
ti = cbind(t1, t2[,2:3])
colnames(ti) = c("ochr", "beg1", "end1", "beg2", "end2")
nrow(ti)
ti = ti[ti$beg2 - ti$end1 + 1 < 0.1*(ti$end1-ti$beg1+1), ]
nrow(ti)

ti = merge(ti, tm[,c(1,3)], by = 'ochr')
ti = cbind(ti, len1 = ti$end1 - ti$beg1 + 1, len2 = ti$end2 - ti$beg2 + 1)
gri1 = with(ti, GRanges(seqnames = nchr, ranges = IRanges(beg1, end = end1)))
gri2 = with(ti, GRanges(seqnames = nchr, ranges = IRanges(beg2, end = end2)))

ocntn1 = intersect_count(gri1, grgn)
ocntn2 = intersect_count(gri2, grgn)
ocntc1 = intersect_count(gri1, grgc)
ocntc2 = intersect_count(gri2, grgc)

ti = cbind(ti, onbs1 = ocntn1, ocrp1 = ocntc1, onbs2 = ocntn2, ocrp2 = ocntc2)

fv = sprintf("%s/%s_HM101/31_sv/01.stb", Sys.getenv("misc3"), org)
tv = read.table(fv, header = T, sep = "\t", as.is = T)
tvs = tv[tv$qlen > 50,]
grv = with(tvs, GRanges(seqnames = qchr, ranges = IRanges(qbeg, end = qend)))

olen1 = intersect_basepair(gri1, grv)
olen2 = intersect_basepair(gri2, grv)

to = cbind(ti, olen1 = olen1, olen2 = olen2)
to = to[(to$olen1+to$olen2) / (to$len1+to$len2) >= 0.3,]
to = to[order(to$len1, decreasing = T), ]
nrow(to)

tos = to[to$onbs1 > 0 | to$ocrp1 > 0,]
tos = tos[order(tos$ocrp1 + tos$onbs1, decreasing = T),]
nrow(tos)
sum(tos$len1 >= 3000)

## check ovlp with SV and synteny blocks
fg = file.path(dirg, "51.gtb")
tg = read.table(fg, header = T, sep = "\t", as.is = T)
tg = tg[,c(1,3:6,16,17)]

grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
grgn = grg[tg$cat2 %in% c("TIR-NBS-LRR", "CC-NBS-LRR")]
grgc = grg[tg$cat2 %in% c("CRP0000-1030", "NCR", "CRP1600-6250")]

dg = tg[tg$cat2 %in% c("CRP0000-1030", "NCR", "CRP1600-6250"),]

fv = sprintf("%s/%s_HM101/31_sv/01.stb", Sys.getenv("misc3"), org)
tv = read.table(fv, header = T, sep = "\t", as.is = T)
tvs = tv[tv$qlen > 50,]
grv = with(tvs, GRanges(seqnames = qchr, ranges = IRanges(qbeg, end = qend)))

fy = sprintf("%s/%s_HM101/23_blat/31.9/gax", Sys.getenv("misc3"), org)
ty = read.table(fy, header = F, sep = "\t", as.is = T)
colnames(ty) = c("tchr", "tbeg", "tend", "tsrd", "qchr", "qbeg", "qend", "qsrd", "cid", "lev")
tys = ty[ty$lev <= 2,]
gry = with(tys, GRanges(seqnames = qchr, ranges = IRanges(qbeg, end = qend)))

len = dg$end - dg$beg + 1
volen = intersect_basepair(grgc, grv)
yolen = intersect_basepair(grgc, gry)

dg2 = cbind(dg, vop = volen / len, yop = yolen / len)
dg3 = dg2[dg2$yop < 0.5 | dg2$vop > 0.5, c(-5,-6)]
nrow(dg3)

fo = sprintf("%s/%s.tbl", dirw, org)
write.table(dg3, file = fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')



#####
source("comp.plot.fun.R")
tracks = c('tgene', 'taxis', 'tgap', 'link', 'qgap', 'qaxis', 'qgene', 'qpacbio')
tracks = c('tgene', 'taxis', 'tgap', 'link', 'qgap', 'qaxis', 'qgene')

fl = file.path(dirw, 'loci.xlsx')
tl = read.xlsx(fl, sheetIndex = 1, header = T)

i = 14
tls = tl[tl$idx == i,]
gro =  with(tls, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

qnames = strsplit(as.character(tls$qnames), split = ' ')[[1]]

cfgs = get_genome_cfgs(c(tname, qnames))
dats = prep_plot_data(gro, cfgs, tname, qnames, tracks)
res = comp.plot(dats, tname, qnames, tracks, draw.title = T)

fn = sprintf("%s/figs/%04d.pdf", dirw, i)
pdf(file = fn, width = 10, height = res$ht/72, bg = 'transparent')
grid.newpage()
grid.draw(res$grobs)
dev.off()