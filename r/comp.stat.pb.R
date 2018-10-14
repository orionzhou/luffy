require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(gtable)
require(RColorBrewer)
source('Location.R')
source('comp.fun.R')

diro = file.path(Sys.getenv("misc2"), "pb.stat")

qnames = c("HM340", "HM340.PB", "HM340.PBBN", "HM340.PBDT", "HM340.PBBNDT", "HM340.PBDTBN", "HM340.FN")


tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

tgen = read.table(file.path(tcfg$gdir, "51.gtb"), sep = "\t", header = T, as.is = T)
tg = data.frame(chr = tgen$chr, beg = tgen$beg, end = tgen$end)
grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

trep = read.table(file.path(tcfg$gdir, "12.rm.bed"), sep = "\t", header = F, as.is = T)
tr = data.frame(chr = trep$V1, beg = trep$V2+1, end = trep$V3)
grr = with(tr, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

##### functional annotation stat
stats = list()
statd = list()
nstats = list()
rstats = list()
cstats = list()
tstats = list()

orgs_rnaseq = c("HM034", "HM056", "HM101", "HM340")

fi = file.path(Sys.getenv("misc2"), "genefam", "nbs.info")
nfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("misc2"), "genefam", "crp.info")
cfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("misc2"), "genefam", "rlk.info")
rfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("data"), 'db', 'pfam', 'genefam.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
tfams = ti$dom[ti$fam == 'TE']
fams = unique(ti$fam[order(ti$pri)])

for (qname in c(tname, qnames_15)) {
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  fg = file.path(dir, "51.gtb")
  tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,15:17)]
  ngene = nrow(tg)
  ids_nonte = tg$id[tg$cat2 != "TE"]
  n_nonte = length(ids_nonte)
  
  dtn = table(tg$cat2)
  x = c()
  x[fams] = 0
  x[names(dtn)] = dtn
  statd[[qname]] = matrix(c(x, ngene),
    nrow = 1, dimnames = list(NULL, c(fams, "Total Genes")))
  
  tf = file.path(dir, "51.fas")
  y = read.fasta(tf, as.string = TRUE, seqtype = "AA")
  dl = ldply(y, nchar)
  lens = dl$V1[dl$.id %in% ids_nonte]
  mean_prot = mean(lens)
  med_prot = median(lens)

  tt = tg[tg$cat2 == 'TE',]
  dtn = table(tt$cat3)
  dtn = dtn[names(dtn) != '']
  x = c()
  x[tfams] = 0
  x[names(dtn)] = dtn
  tstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, tfams))
  n_te = nrow(tt)
  
  tn = tg[tg$cat2 %in% c('CC-NBS-LRR', 'TIR-NBS-LRR'),]
  dtn = table(tn$cat3)
  x = c()
  x[nfams] = 0
  x[names(dtn)] = dtn
  nstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, nfams))
  n_nbs = nrow(tn)
  
  tr = tg[tg$cat2 %in% c('RLK','LRR-RLK'),]
  dtn = table(tr$cat3)
  x = c()
  x[rfams] = 0
  x[names(dtn)] = dtn
  rstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, rfams))
  n_rlk = nrow(tr)
  
  tc = tg[tg$cat2 %in% c('CRP0000-1030','NCR','CRP1600-6250'),]
  dtn = table(tc$cat3)
  x = c()
  x[cfams] = 0
  x[names(dtn)] = dtn
  cstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, cfams))
  n_ncr = sum(tc$cat2 == 'NCR')
  
  tf = tg[tg$cat2 %in% c('F-box'),]
  n_fbx = nrow(tf)
  
  ids = ids_nonte
  ids_rna = c(); ids_hom = c()
  if(qname %in% orgs_rnaseq) {
    org = qname
    ff = file.path(Sys.getenv("misc2"), "rnaseq/mt/31_cufflinks", org, 'isoforms.fpkm_tracking')
    tf = read.table(ff, header = T, sep = "\t", as.is = T)
    ids_rna = tf$tracking_id[tf$tracking_id %in% ids & tf$FPKM > 0]
    p_rna = sprintf("%.01f", length(ids_rna) / length(ids) * 100)
  } else {
    p_rna = '-'
  }
  if(qname != tname) {
    ft = sprintf("%s/%s_HM101/51_ortho/31.ortho.tbl", Sys.getenv("misc3"), qname)
    tt = read.table(ft, header = T, sep = "\t", as.is = T)
    ids_hom = tt$qid[tt$qid %in% ids]
    p_hom = sprintf("%.01f", length(ids_hom) / length(ids) * 100)
  } else {
    p_hom = '-'
  }
  p_sup = sprintf("%.01f", length(unique(c(ids_rna, ids_hom))) / length(ids) * 100)
  
  stats[[qname]] = matrix(c(ngene, n_te, n_nonte, n_nbs, n_fbx, n_rlk, n_ncr, med_prot, p_rna, p_hom, p_sup), nrow = 1, dimnames = list(NULL, c("# Total Genes", 'TE', 'non-TE','NBS-LRR','F-box','LRR-RLK','NCR',"Median Prot Length", 'RNA-seq (%)', 'Homology (%)', 'RNA-seq + Homology (%)')))
}
ds = do.call(rbind.data.frame, stats)
#for (i in 1:ncol(do)) { do[,i] = format(do[,i], big.mark = ",") }
fo = file.path(diro, "03_annotation.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)

ds = do.call(rbind.data.frame, statd)
fo = file.path(diro, "04.annotation.detail.tbl")
write.table(t(ds), fo, sep = "\t", row.names = T, col.names = T, quote = F)

ds = do.call(rbind.data.frame, nstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "04.nbs.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

ds = do.call(rbind.data.frame, rstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "04.rlk.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

ds = do.call(rbind.data.frame, cstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "04.crp.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

ds = do.call(rbind.data.frame, tstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "04.te.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

##### comparative / novel-seq stats
stats = list()
for (qname in qnames) {
  #add repeatmasker stats
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  flen = file.path(dir, "15.sizes")
  tlen = read.table(flen, sep = "\t", header = F, as.is = T)
  total_len = sum(tlen$V2)
  
  fgap = file.path(dir, "16.gap.bed")
  if(file.info(fgap)$size == 0) {
  	tgap = data.frame(V1=c(),V2=c(),V3=c())
  } else {
  	tgap = read.table(fgap, sep = "\t", header = F, as.is = T)
  }
  total_gap = sum(tgap$V3 - tgap$V2)
  total_bases = total_len - total_gap
  
  fcds = file.path(dir, "51.tbl")
  tcds = read.table(fcds, sep = "\t", header = F, as.is = T)[,1:6]
  tcds = tcds[tcds$V6 == 'cds',]
  grc = GRanges(seqnames = tcds$V1, ranges = IRanges(tcds$V2, end = tcds$V3))
  grc = reduce(grc)
  
  frep = file.path(dir, "12.rm.bed")
  trep = read.table(frep, sep = "\t", header = F, as.is = T)
  brep = sum(trep$V3 - trep$V2)
  pct_rep = brep / total_bases * 100

  dir = sprintf("%s/%s_%s/23_blat", Sys.getenv("misc3"), qname, tname)
  fi = file.path(dir, '41.5/gax')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:3,10)]
  colnames(ti) = c('tid', 'tbeg', 'tend', 'lev')
  aligned = sum(ti$tend - ti$tbeg + 1)
#  pct_aligned = aligned / total_bases * 100

  fi = file.path(dir, '31.9/gax')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:3,10)]
  colnames(ti) = c('tid', 'tbeg', 'tend', 'lev')
  ti = ti[ti$lev <= 2,]
  synteny = sum(ti$tend - ti$tbeg + 1)
#  pct_synteny = synteny / total_bases * 100

  dir = sprintf("%s/%s_%s/41_novseq", Sys.getenv("misc3"), qname, tname)
  fi = file.path(dir, '21.bed')
  ti = read.table(fi, sep = "\t", header = F, as.is = T)
  bnov = sum(ti$V3 - ti$V2)
  pnov = bnov / total_bases * 100
  
  grn = GRanges(seqnames = ti$V1, ranges = IRanges(ti$V2+1, end = ti$V3))
  grn = reduce(grn)
  bcds = sum(width(intersect(grc, grn)))
  pcds = bcds / bnov * 100

  total_bases = format(total_bases, big.mark = ",")
  brep = format(brep, big.mark = ",")
#  pct_rep = sprintf("%.01f%%", pct_rep)
  aligned = format(aligned, big.mark = ",")
#  pct_aligned = sprintf("%.01f%%", pct_aligned)
  synteny = format(synteny, big.mark = ",")
#  pct_synteny = sprintf("%.01f%%", pct_synteny)
  bnov = format(bnov, big.mark = ",")
  pnov = sprintf("%.01f%%", pnov)
  bcds = format(bcds, big.mark = ",")
  pcds = sprintf("%.01f%%", pcds)
  stats[[qname]] = matrix(c(total_bases, brep, aligned, synteny, bnov, pnov, bcds, pcds), 
    nrow = 1, dimnames = list(NULL, c("Total Bases", "Repetitive", "Alignable to HM101", "Bases in Synteny Blocks", "Novel Sequences", "", "Novel Coding Seq", "")))
}
ds = do.call(rbind.data.frame, stats)
fo = file.path(diro, "11_comp_stat.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)

### synteny heatmap
chrn = 1:nrow(tt)
names(chrn) = rev(tt$chr)

qname = "HM340.FN"
qname = "HM034"

dirq = sprintf("%s/%s", Sys.getenv("genome"), qname)
fl = file.path(dirq, "15.sizes")
tl = read.table(fl, header = F, sep = "\t", as.is = T)
colnames(tl) = c("chr", "len")

dirw = sprintf("%s/%s_HM101/23_blat", Sys.getenv("misc3"), qname)
fi = file.path(dirw, "31.9.gal")
ti = read.table(fi, header = T, sep = "\t", as.is = T)[,1:15,]

tt2 = cbind(tt, chrn = chrn[as.character(tt$chr)])
ti2 = cbind(ti, chrn = chrn[as.character(ti$tId)])
tg2 = cbind(tg, chrn = chrn[as.character(tg$chr)])

ti3 = cbind(ti2, lev2 = ti2$lev, y=ti2$chrn)
ti3$lev2[ti3$lev == 1] = '1'
ti3$lev2[ti3$lev == 2] = '2'
ti3$lev2[ti3$lev == 3] = '3'
ti3$lev2[ti3$lev >= 4] = '4'
ti3$y[ti3$lev == 1] = ti3$y[ti3$lev == 1] + 0.15
ti3$y[ti3$lev == 2] = ti3$y[ti3$lev == 2] + 0.05
ti3$y[ti3$lev == 3] = ti3$y[ti3$lev == 3] - 0.05
ti3$y[ti3$lev >= 4] = ti3$y[ti3$lev >= 4] - 0.15
breaks = c('1','2','3', '4')
cols = brewer.pal(4, "Set2")

chr_coord <- function(l) {
     sprintf("%dMb", floor(l/1000000))
}
pc <- ggplot() +
  geom_rect(data = tt2, aes(xmin=beg, xmax=end, ymin=chrn-0.3, ymax=chrn+0.3), fill=NA, color='black') +
  geom_rect(data = ti3, aes(xmin=tBeg, xmax=tEnd, ymin=y-0.05, ymax=y+0.05, fill=lev2), color=NA) +
  geom_rect(data = tg2, aes(xmin=beg, xmax=end, ymin=chrn+0.2, ymax=chrn+0.3), fill='royalblue', color=NA) +
#  geom_point(data = ti2, aes(x=tBeg, y=chrn+0.4, col=qId, shape=qId), size=2) +
#  geom_point(data = tm2, aes(x=tBeg, y=chrn+0.4, col='Mt3.5 Cen', shape='Mt3.5 #Cen'), size=2) +
#  geom_point(data = tc2, aes(x=pos, y=chrn-0.35, col='Predicted Mt4.0 Cen', shape='Predicted Mt4.0 Cen'), size=2) +
  theme_bw() + 
  scale_x_continuous(name = '', expand = c(0.01, 0), labels=chr_coord) + 
  scale_y_continuous(name = '', expand = c(0, 0), breaks=chrn, labels=names(chrn), limits=c(0.5,9.7)) +
  scale_fill_manual(name='', breaks = breaks, values=cols) +
#  scale_color_manual(name='', values=c('red','blue','green','purple','black')) +
#  scale_shape_manual(name='', values=c(6,6,6,6,17)) +
#  theme(legend.position = c(0.9,0.45), legend.direction = "vertical", legend.title = element_text(size = 8), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size=8), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "cm")) +
#  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))

fo = sprintf("%s/51.syn.%s.pdf", diro, qname)
ggsave(pc, filename = fo, width = 10, height = 6)


## look for TE and lev2 syn-block
tg = data.frame(chr = tgen$chr, beg = tgen$beg, end = tgen$end)[tgen$cat2=='TE',]
grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

tx = data.frame(chr = ti$tId, beg = ti$tBeg, end = ti$tEnd)[ti$lev==2,]
grx = with(tx, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

bp = intersect_basepair(grx, union(grr,grg))
nrow(tx)
sum(bp/width(grx) > 0)
sum(bp/width(grx) > 0.5)

for (qname in c("HM034", "HM340", "HM340.PB", "HM340.PBBN", "HM340.PBDT", "HM340.PBBNDT", "HM340.PBDTBN")) {
	xxx(qname)
}

#### visualize with centromere locations
qname = "HM340.FN"
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

tt = tt[1:49,]
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
breaks1 = seq(vmin, vmax, length.out = length(cols1)+1)
gapcol = cut(gap, breaks1, include.lowest = T)

vmin = min(rep); vmax = max(rep)
breaks2 = seq(vmin, vmax, length.out = length(cols2)+1)
repcol = cut(rep, breaks1, include.lowest = T)

tw = cbind(tw, gapcol=gapcol, repcol=repcol)

## cent
dirw = file.path(Sys.getenv("misc2"), "centromere")
fc = file.path(dirw, "51.HM340.gal")
tc = read.table(fc, header = T, sep = "\t", as.is = T)
tc = tc[tc$score > 100,c(2:11,13:14,18:19)]
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
  scale_y_continuous(name = '', expand = c(0, 0), breaks=chrn, labels=names(chrn), limits=c(0.5,49.7)) +
  scale_fill_manual(name='', breaks=breaks1, values=cols1, guide=guide_legend(label.position='none')) +
  scale_color_manual(name='', values=c('red','blue','black','purple','red')) +
  scale_shape_manual(name='', values=c(6,6,6,6,17)) +
  theme(legend.position = c(0.8,0.8), legend.direction = "vertical", legend.title = element_text(size = 8), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size=8), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "cm")) +
#  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))

fo = file.path(dirw, "92.pdf")
ggsave(pc, filename = fo, width = 8, height = 8)


