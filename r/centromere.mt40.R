require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(gtable)
require(grid)
require(RColorBrewer)
source('Location.R')
source('comp.fun.R')

dirw = file.path(Sys.getenv("misc2"), "centromere")

tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2, stringsAsFactors=F)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3, stringsAsFactors=F)
tp = cbind(tp, len = tp$end-tp$beg+1)

### find largest gaps - not working
xf <- function(df) {
  df = cbind(df, len = df[,'end']-df[,'beg']+1)
  idxs = which(df[,'len']==max(df[,'len']))
  data.frame(chr=df[idxs,'chr'], beg=df[idxs,'beg'], end=df[idxs,'end'], len=df[idxs,'len'], stringsAsFactors=F)
}
tl = ddply(tp, .(chr), xf)
tl = tl[tl$chr != 'chrU',]
fl = file.path(dirw, "31.lgap.mt40.tbl")
write.table(tl, fl, sep = "\t", row.names = F, col.names = T, quote = F)

### using Mt3.5 centromeres
fi = file.path(Sys.getenv("genome"), "Mtruncatula_3.5/sequence/07_cen.tbl")
ti = read.table(fi, header = F, sep = "\t", as.is = T)
dc1 = within(ti, {V3=V2-1; V2 = V2-1000;})
dc2 = within(ti, {V2=V3+1; V3 = V3+1000;})
dm = rbind(dc1, dc2)
fm = file.path(dirw, "61.cen.mt35.bed")
write.table(dm, fm, sep = "\t", row.names = F, col.names = F, quote = F)
# seqret.pl -d $genome/Mtruncatula_3.5/11_genome.fa -b 61.cen.mt35.bed -o 61.fas
# blat $genome/HM101/11_genome.fas 61.fas 62.psl
# psl2gal.pl -i 62.psl -o 62.gal

fi = file.path(dirw, "62.gal")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
qChr = sapply(strsplit(ti$qId, "-"), "[", 1) 
tm = ti[ti$mat > 950 & ti$tId == qChr,]

### prepare sliding window of gap stat
x = tt$end
names(x) = tt$chr
gr = tileGenome(x, tilewidth = 100000, cut.last.tile.in.chrom = T)

tw = data.frame(chr = seqnames(gr), beg = start(gr), end = end(gr), 
  len = width(gr), stringsAsFactors = F)
bp_gap = intersect_basepair(gr, grp)
gap = bp_gap / tw$len

cols = brewer.pal(9, "Greys")

vmin = min(gap); vmax = max(gap)
breaks = seq(vmin, vmax, length.out = length(cols)+1)
gapcol = cut(gap, breaks, include.lowest = T)
tw = cbind(tw, gap=gap, gap=gapcol)


chrn = 1:nrow(tt)
names(chrn) = rev(tt$chr)

fi = file.path(dirw, "11.gal")
ti = read.table(fi, header = T, sep = "\t", as.is = T)

fc = file.path(dirw, "cen.mt40.tsv")
tc = read.table(fc, header = F, sep = "\t", as.is = T)
colnames(tc) = c('chr','beg','end')
tc = cbind(tc, pos = (tc$end+tc$beg)/2)

tt2 = cbind(tt, chrn = chrn[as.character(tt$chr)])
tw2 = cbind(tw, chrn = chrn[as.character(tw$chr)])
ti2 = cbind(ti, chrn = chrn[as.character(ti$tId)])
tm2 = cbind(tm, chrn = chrn[as.character(tm$tId)])
tc2 = cbind(tc, chrn = chrn[as.character(tc$chr)])

pc <- ggplot() +
  geom_rect(data = tw2, aes(xmin=beg, xmax=end, ymin=chrn-0.25, ymax=chrn+0.25, fill=gapcol), linetype=0) +
  geom_rect(data = tt2, aes(xmin=beg, xmax=end, ymin=chrn-0.25, ymax=chrn+0.25), fill=NA, color='deeppink4') +
  geom_point(data = ti2, aes(x=tBeg, y=chrn+0.4, col=qId, shape=qId), size=2) +
  geom_point(data = tm2, aes(x=tBeg, y=chrn+0.4, col='Mt3.5 Cen', shape='Mt3.5 Cen'), size=2) +
  geom_point(data = tc2, aes(x=pos, y=chrn-0.35, col='Predicted Mt4.0 Cen', shape='Predicted Mt4.0 Cen'), size=2) +
  theme_bw() + 
  scale_x_continuous(name = '', expand = c(0.01, 0)) + 
  scale_y_continuous(name = '', expand = c(0, 0), breaks=chrn, labels=names(chrn), limits=c(0.5,9.7)) +
  scale_fill_manual(name='', breaks=breaks, values=cols, guide=guide_legend(label.position='none')) +
  scale_color_manual(name='', values=c('red','blue','green','purple','black')) +
  scale_shape_manual(name='', values=c(6,6,6,6,17)) +
  theme(legend.position = c(0.9,0.45), legend.direction = "vertical", legend.title = element_text(size = 8), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size=8), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "cm")) +
#  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))

fo = file.path(dirw, "91.pdf")
ggsave(pc, filename = fo, width = 10, height = 6)


## remix
qname = "HM101"
dirg = sprintf("%s/%s", Sys.getenv("genome"), qname)
flen = file.path(dirg, "15.bed")
tt = read.table(flen, sep = "\t", header = F, as.is = T)
colnames(tt) = c("chr", "beg", "end")
tt$beg = tt$beg + 1

fgap = file.path(dirg, "16.gap.bed")
tgap = read.table(fgap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

fgen = file.path(dirg, "51.gtb")
tgen = read.table(fgen, sep = "\t", header = T, as.is = T)
tg = data.frame(chr = tgen$chr, beg = tgen$beg, end = tgen$end)
grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

frep = file.path(dirg, "12.rm.bed")
trep = read.table(frep, sep = "\t", header = F, as.is = T)
tr = data.frame(chr = trep$V1, beg = trep$V2+1, end = trep$V3)
grr = with(tr, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

tt = tt[1:8,]
chrs = tt$chr

chrn = 1:nrow(tt)
names(chrn) = rev(tt$chr)

## gap & repeat stat
x = tt$end
names(x) = tt$chr
gr = tileGenome(x, tilewidth = 100000, cut.last.tile.in.chrom = T)

tw = data.frame(chr = seqnames(gr), beg = start(gr), end = end(gr), 
  len = width(gr), stringsAsFactors = F)
bp_gap = intersect_basepair(gr, grp)
bp_rep = intersect_basepair(gr, grr)
rep = bp_rep / tw$len
gap = bp_gap / tw$len

cols1 = brewer.pal(9, "Greys")
cols2 = brewer.pal(9, "OrRd")

vmin = min(gap); vmax = max(gap)
breaks1 = seq(0, 1, length.out = length(cols1)+1)
gapcol = cut(gap, breaks1, include.lowest = T)

vmin = min(rep); vmax = max(rep)
breaks2 = seq(0, 1, length.out = length(cols2)+1)
repcol = cut(rep, breaks1, include.lowest = T)

tw = cbind(tw, gapcol=gapcol, repcol=repcol)

## cent
dirw = file.path(Sys.getenv("misc2"), "centromere")
fc = file.path(dirw, "11.gal")
tc = read.table(fc, header = T, sep = "\t", as.is = T)
tc = tc[tc$score > 100, c(2:11,13:14,18:19)]
tc2 = ddply(tc, .(tId, qId), summarise, beg = min(tBeg), end=max(tEnd), nrep = length(tId))
tc2 = tc2[tc2$nrep > 5 & tc2$tId %in% chrs,]
tc2 = cbind(tc2, pos = (tc2$beg+tc2$end)/2)

tt2 = cbind(tt, chrn = chrn[as.character(tt$chr)])
tw2 = cbind(tw, chrn = chrn[as.character(tw$chr)])
tc2 = cbind(tc2, chrn = chrn[as.character(tc2$tId)])

chr_coord <- function(l) {
     sprintf("%dMb", floor(l/1000000))
}
pc <- ggplot() +
  geom_rect(data = tw2, aes(xmin=beg, xmax=end, ymin=chrn, ymax=chrn+0.35, fill=gapcol), linetype=0) +
  geom_rect(data = tw2, aes(xmin=beg, xmax=end, ymin=chrn-0.35, ymax=chrn, fill=repcol), linetype=0) +
  geom_rect(data = tt2, aes(xmin=beg, xmax=end, ymin=chrn-0.35, ymax=chrn+0.35), fill=NA, color='black') +
  geom_segment(data = tt2, aes(x=beg, xend=end, y=chrn, yend=chrn), color='black') +
  geom_point(data = tc2, aes(x=pos, y=chrn+0.5, col='MtR3 Centromeric Repeat', shape='MtR3 Centromeric Repeat'), size=1.5) +
  theme_bw() + 
  scale_x_continuous(name = '', expand = c(0.01, 0), labels = chr_coord) + 
  scale_y_continuous(name = '', expand = c(0, 0), breaks=chrn, labels=names(chrn), limits=c(0.5,8.7)) +
  scale_fill_manual(name='', breaks=breaks1, values=cols1, guide=guide_legend(label.position='none')) +
  scale_color_manual(name='', values=c('red','blue','black','purple','red')) +
  scale_shape_manual(name='', values=c(6,6,6,6,17)) +
  theme(legend.position = c(0.8,0.3), legend.direction = "vertical", legend.title = element_text(size = 8), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size=8), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "cm")) +
#  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))

fo = file.path(dirw, "90.pdf")
ggsave(pc, filename = fo, width = 5, height = 5)
