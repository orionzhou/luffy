require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(gtable)
require(grid)
require(gridExtra)
require(RColorBrewer)
#require(Gviz)
source('Location.R')
source('comp.win.fun.R')

dirw = file.path(Sys.getenv("misc3"), "comp.stat")

##### read in data tracks
tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

### gene fam stats
tg = read.table(tcfg$gene, sep = "\t", header = F, as.is = T)[,1:7]
colnames(tg) = c("chr", 'beg', 'end', 'srd', 'id', 'type', 'cat')
tg = tg[tg$type == 'cds',]

### covered bases
grl = list()
for (qname in qnames_15) {
  dirc = sprintf("%s/%s_%s", Sys.getenv("misc3"), qname, tname) 
  fa = file.path(dirc, '23_blat/31.9/gax')
  ta = read.table(fa, header = T, sep = "\t", as.is = T)
  colnames(ta) = c('tchr','tbeg','tend','tsrd','qchr','qbeg','qend','qsrd','cid','lev')
  gra = with(ta[ta$lev<=20,], GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend)))
  grl[[qname]] = reduce(gra)
}

### ingroup coverage & theta-pi
fc = file.path(Sys.getenv("misc3"), "comp.vnt", "81.cvg.tbl")
tc = read.table(fc, header = F, sep = "\t", as.is = T)
grc = with(tc, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

fn = file.path(Sys.getenv("misc3"), "comp.vnt", "25.stat.tbl")
tn = read.table(fn, header = T, sep = "\t", as.is = T)
grn = with(tn, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))

fv = file.path(Sys.getenv("misc3"), "comp.vnt", "52.stat.tbl")
tv = read.table(fv, header = T, sep = "\t", as.is = T)
tv = cbind(tv, size = (tv$rsize + tv$asize - 2))
tvs = tv[tv$size < 50,]
tvl = tv[tv$size >= 50,]
grvs = with(tvs, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))
grvl = with(tvl, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))

##### create sliding window table using 10kb sliding windows
chr = 'chr5'
beg = 7
title = sprintf("%s:%02dMb", chr, beg)

res = prepare_data(chr, beg, grt, grp, tg, grl, grc, grn, grvs, grvl)
splots = sub_plots(chr, beg, res$tw, res$dg, res$dy, res$ds)

p_syn = splots$syn + theme(plot.margin = unit(c(0.1,0.1,0.1,1), "lines"))
p_sta = splots$sta
gr_syn = ggplotGrob(p_syn)
gr_sta = ggplotGrob(p_sta)
gr_sta$widths = gr_syn$widths
#gt1 = gtable_add_cols(gt1, unit(0, "mm"))
gs = list(gr_syn, gr_sta)
heis = c(3, 2)

g <- gtable_matrix(name='demo', grobs = matrix(gs, nrow = length(gs)), widths = 1, heights = heis)
fo = sprintf("%s/33.wins/%s.pdf", dirw, title)
pdf(file = fo, width = 3, height = 6, bg = 'transparent')
grid.newpage()
grid.draw(g)
dev.off()

##### create multiple sliding-window plot
source('comp.win.fun.R')
chrs = c('chr2', 'chr2', 'chr4', 'chr5', 'chr7')
begs = c(16, 30, 5, 6, 28)
chrs = c('chr7', 'chr2', 'chr6', 'chr2', 'chr8')
begs = c(28, 16, 5, 30, 25)

lres = list()
for (i in 1:length(chrs)) {
  res = prepare_data(chrs[i], begs[i], grt, grp, tg, grl, grc, grn, grvs, grvl)
  lres[[i]] = res
}

gts = list()
for (i in 1:length(chrs)) {
  res = lres[[i]]
  splots = sub_plots(chrs[i], begs[i], res$tw, res$dg, res$dy, res$ds)
  p_syn = splots$syn + theme(plot.margin = unit(c(0.5,0.1,0,0.5), "lines"))
  p_syn = p_syn + theme(plot.title = element_blank())
  p_sta = splots$sta + theme(plot.margin = unit(c(0.2,0.1,0.1,0.5), "lines"))
  
  if(i == 1) {
    p_syn = p_syn + theme(plot.margin = unit(c(0.5,0.1,0,1), "lines"))
    p_sta = p_sta + theme(plot.margin = unit(c(0.2,0.1,0.1,1), "lines"))
  } else {
    p_syn = p_syn + theme(axis.text.y = element_blank())
    p_sta = p_sta + theme(axis.text.y = element_blank())
  }
  gr_syn = ggplotGrob(p_syn)
  gr_sta = ggplotGrob(p_sta)
  gr_sta$widths = gr_syn$widths
  gs = list(gr_syn, gr_sta)
  heis = c(3, 2)
  g <- gtable_matrix(name='demo', grobs = matrix(gs, nrow = length(gs)), widths = 1, heights = heis)
  g_rect = rectGrob(gp = gpar(col='black', fill=NA, lwd=2))
  g_titl = textGrob(LETTERS[i], x=unit(0.5,'npc'), y=unit(1,'npc')-unit(0.1,'lines'), vjust=1)
  g = gtable_add_grob(g, g_titl, t=1, l=1, b=2, r=1)
  gts[[i]] = g
}

bog <- rectGrob(gp = gpar(col='black', fill=NA, lwd=2))
gf <- gtable_matrix(name='demo', grobs = matrix(gts, ncol = length(gts)), widths = c(3.9,rep(3,length(gts)-1)), heights = 1)
fo = sprintf("%s/33.wins/fig2.pdf", dirw)
pdf(file = fo, width = 7, height = 6, bg = 'transparent')
grid.newpage()
grid.draw(gf)
dev.off()

