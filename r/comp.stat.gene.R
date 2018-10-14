require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(gtable)
#require(Gviz)
source('Location.R')
source('comp.fun.R')

diro = file.path(Sys.getenv("misc3"), "comp.stat")

qname = 'HM004'
cfg = cfgs[[qname]]

tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

### enrichment of TEs in non-syntenic regions
fx = file.path(cfg$cdir, "31.9/gax")
tx = read.table(fx, sep = "\t", header = T, as.is = T)
colnames(tx) = c("tchr", "tbeg", "tend", "tsrd", "qchr", "qbeg", "qend", "qsrd", "cid", "lev")
tx = tx[tx$lev == 1,]

fg = file.path(tcfg$gdir, "51.gtb")
dg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,16)]
grg = with(dg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

gs = with(tx, GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend)))

olens = intersect_basepair(grg, gs)
to = cbind(dg, ovlp = olens/(dg$end - dg$beg + 1))
ddply(to, .(cat2), summarise, pct_syn = sum(ovlp > 0.5) / length(ovlp))


fg = file.path(cfg$gdir, "51.gtb")
dg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,16)]
grg = with(dg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

gs = with(tx, GRanges(seqnames = qchr, ranges = IRanges(qbeg, end = qend)))

olens = intersect_basepair(grg, gs)
to = cbind(dg, ovlp = olens/(dg$end - dg$beg + 1))
ddply(to, .(cat2), summarise, pct_syn = sum(ovlp > 0.5) / length(ovlp))

