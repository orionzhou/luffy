require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(rtracklayer)
require(grid)
source("Align.R")
source('Location.R')
source('comp.fun.R')

dirw = file.path(Sys.getenv("misc3"), "comp.syn")

### use liftOver to identify syn-ortho
fg = file.path(Sys.getenv("genome"), tname, "51.tbl")
tg = read.table(fg, header = F, sep = "\t", as.is = T)
colnames(tg) = c('chr','beg','end','srd','id','type','fam')
tg = tg[tg$type == 'cds', c('chr','beg','end','id','fam')]
grt = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gb = group_by(tg, id)
tgl = summarise(gb, len = sum(end - beg + 1))

qnames = qnames_15
for (qname in qnames) {

fg = file.path(Sys.getenv("genome"), qname, "51.tbl")
qg = read.table(fg, header = F, sep = "\t", as.is = T)
colnames(qg) = c('chr','beg','end','srd','id','type','fam')
qg = qg[qg$type == 'cds', c('chr','beg','end','id','fam')]
grq = with(qg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gb = group_by(qg, id)
qgl = summarise(gb, len = sum(end - beg + 1))

fchain = sprintf("%s/%s_HM101/23_blat/31.9.chain", Sys.getenv("misc3"), qname)
chain = import.chain(fchain)

grt = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
y = liftOver(grt, chain)
dm = as(y, 'data.frame')[,c(1,3:6)]
gr1 = with(dm, GRanges(seqnames = seqnames, ranges = IRanges(start, end = end)))

dm2 = cbind(idx1 = dm$group, idx2 = 1:nrow(dm), dm[,2:5])
colnames(dm2)[3:6] = c("cchr", "cbeg", "cend", "clen")

x = intersect_gene(gr1, grq, qg$id)
colnames(x) = c('idx2', 'qid', 'olen')
dx2 = merge(x, qgl, by.x='qid', by.y='id')
colnames(dx2)[4] = 'qlen'

dm3 = merge(dm2, dx2, by='idx2', all.x=T)

dg = cbind(idx1 = 1:nrow(tg), tg[,1:4])
dg = merge(dg, tgl, by='id')
colnames(dg)[1:6] = c("tid", "idx1", "tchr", "tbeg", "tend", 'tlen')
dg2 = merge(dg, dm3, by = 'idx1', all.x = T)

dg3 = dg2[!is.na(dg2$qid),]
gb = group_by(dg3, tid, qid)
dg4 = summarise(gb, tlen=tlen[1], qlen=qlen[1], olen=sum(olen))

dg5 = dg4[dg4$olen / min(dg4$tlen, dg4$qlen) >= 0.4 & dg4$olen >= 10,]
gb = group_by(dg5, tid)
dg6 = summarise(gb, qid=qid[which(olen==max(olen))[1]], tlen=tlen[1], qlen=qlen[which(olen==max(olen))[1]], olen=olen[which(olen==max(olen))[1]])

fo = sprintf("%s/01.syn.pair/%s.tsv", dirw, qname)
write.table(dg6, fo, sep = "\t", row.names = F, col.names = T, quote = F)

}

### calculate similarity by pairwise alignment
cl = makeCluster(detectCores())

cluster_fun <- function() {
    require(Biostrings)
    require(plyr)
}
clusterCall(cl, cluster_fun)

f_tfas = file.path(Sys.getenv("genome"), tname, "51.fas")
tfas <- read.fasta(f_tfas, seqtype = "AA", as.string = T, set.attributes = F)
tids = names(tfas)

qnames = qnames_15
for (qname in qnames) {
  f_qfas = file.path(Sys.getenv("genome"), qname, "51.fas")
  qfas <- read.fasta(f_qfas, seqtype = "AA", as.string = T, set.attributes = F)
  qids = names(qfas)

  fi = sprintf("%s/01.syn.pair/%s.tsv", dirw, qname)
  ti = read.table(fi, sep = "\t", header = T, as.is = T)[,1:5]
  ti = ti[ti$qid != '' & ti$qid %in% qids & ti$tid %in% tids,]

  tm = cbind(ti, tseq = as.character(tfas[ti$tid]), 
    qseq = as.character(qfas[ti$qid]), stringsAsFactors = F)

  ptm <- proc.time()
  y = parApply(cl, tm, 1, aa_pw_dist)

  cat(qname, proc.time() - ptm, "\n")
  
  to = cbind(ti, t(y))
  to$qlen = as.integer(to$qlen / 3)
  to$tlen = as.integer(to$tlen / 3)
  fo = sprintf("%s/05.score/%s.tsv", dirw, qname)
  write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
}

stopCluster(cl)

#### create syn-ortho ID and score matrix
fg = file.path(tcfg$gdir, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:5,16)]

tg = tg[order(tg$chr, tg$beg, tg$end),]
tg = cbind(tg, idx = 1:nrow(tg))
colnames(tg)[6] = 'fam'
tg = rename_genefam(tg)

td = tg[,c(1,7)]
colnames(td)[1] = tname
tc = tg[,c(1,7)]
colnames(tc)[1] = tname

for (qname in qnames_15) {
  fs = sprintf("%s/05.score/%s.tsv", dirw, qname)
  tt = read.table(fs, sep = "\t", header = T, as.is = T)
  tt = cbind(tt, cvg = 1 - (tt$qgap + tt$tgap) / tt$len, sim = tt$mat / (tt$mat + tt$mis))
  
  tt2 = tt[tt$cvg >= 0.5 & tt$sim >= 0.4, c('tid','qid')]
  colnames(tt2) = c(tname, qname)
  td = merge(td, tt2, by = tname, all = T)
  
  tt2 = tt[tt$cvg >= 0.5 & tt$sim >= 0.4, c('tid','sim')]
  colnames(tt2) = c(tname, qname)
  tc = merge(tc, tt2, by = tname, all = T)
  
  cat(qname, nrow(tc), nrow(tt2), "\n")
}
td = td[order(td$idx), c(1,3:ncol(td))]
tc = tc[order(tc$idx), c(1,3:ncol(tc))]
stopifnot(sum(is.na(td)) == sum(is.na(tc)))

norgs = apply(td, 1, function(x) sum(!is.na(x[-1])))
table(norgs)

fd = file.path(dirw, "11.ids.tsv")
write.table(td, fd, sep = "\t", row.names = F, col.names = T, quote = F, na = '')
fc = file.path(dirw, "12.score.tsv")
write.table(tc, fc, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

### add insertion/deletion to syn-ortho matrix
fg = file.path(Sys.getenv("genome"), tname, "51.tbl")
tg = read.table(fg, header = F, sep = "\t", as.is = T)
colnames(tg) = c('chr','beg','end','srd','id','type','fam')
tg = tg[tg$type == 'cds', c('chr','beg','end','id','fam')]
grt = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gb = group_by(tg, id)
tgl = summarise(gb, len = sum(end - beg + 1))

qnames = qnames_15
for (qname in qnames) {

fg = file.path(Sys.getenv("genome"), qname, "51.tbl")
qg = read.table(fg, header = F, sep = "\t", as.is = T)
colnames(qg) = c('chr','beg','end','srd','id','type','fam')
qg = qg[qg$type == 'cds', c('chr','beg','end','id','fam')]
grq = with(qg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gb = group_by(qg, id)
qgl = summarise(gb, len = sum(end - beg + 1))

fa = sprintf("%s/%s_HM101/23_blat/31.9.gal", Sys.getenv("misc3"), qname)
ta = read.table(fa, header = T, sep = "\t", as.is = T)[,1:12]

gct = with(ta, GRanges(seqnames = tId, ranges = IRanges(tBeg, end = tEnd)))
gcq = with(ta, GRanges(seqnames = qId, ranges = IRanges(qBeg, end = qEnd)))

x = intersect_gene(gct, grt, tg$id)
x = merge(x, tgl, by.x='qidx', by.y='id')
x = x[x$olen/x$len >= 0.6,]
dct = data.frame(cid = ta$id[x$idx], gid = x$qidx)

x = intersect_gene(gcq, grq, qg$id)
x = merge(x, qgl, by.x='qidx', by.y='id')
x = x[x$olen/x$len >= 0.6,]
dcq = data.frame(cid = ta$id[x$idx], gid = x$qidx)
