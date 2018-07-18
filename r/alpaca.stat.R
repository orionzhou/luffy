require(plyr)
require(dplyr)
require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(xlsx)
require(seqinr)
require(RColorBrewer)
source("Location.R")
source("comp.fun.R")

diro = file.path(Sys.getenv('misc3'), 'alpaca')

qnames = qnames_alpaca_comp

tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

##### functional annotation stat
stats = list()
statd = list()
nstats = list()
rstats = list()
cstats = list()
tstats = list()

fi = file.path(Sys.getenv("misc2"), "genefam", "nbs.info")
nfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("misc2"), "genefam", "crp.info")
cfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("misc2"), "genefam", "rlk.info")
rfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("data"), 'db', 'pfam', 'genefam.tbl')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
tfams = ti$dom[ti$fam == 'TE']
fams = unique(ti$fam[order(ti$pri)])

for (qname in c(tname, qnames)) {
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
  
  tr = tg[tg$cat2 %in% c('LRR-RLK'),]
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
  
  stats[[qname]] = matrix(c(ngene, n_te, n_nonte, n_nbs, n_fbx, n_rlk, n_ncr, med_prot),
    nrow = 1, dimnames = list(NULL, c("# Total Genes", 'TE', 'non-TE', 
    'NBS-LRR','F-box','LRR-RLK','NCR',"Median Prot Length")))
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


### tandem clustering statistics in HM101
fi = sprintf("%s/gene.cluster/11.tandem/%s.tbl", Sys.getenv('misc2'), "HM101")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti = ti[ti$cat2 != 'TE',]
idx_single = which(is.na(ti$clu))
ti$clu[idx_single] = seq(max(ti$clu, na.rm = T)+1, by = 1, length.out = length(idx_single))

brks = c(seq(0.5, 10.5, by = 1), 15.5, Inf)
labs = c(1:10, '11-15', '16+')
labs = factor(labs, levels = labs)

tx = ti
tx$cat2[tx$cat2 %in% c("CC-NBS-LRR", "TIR-NBS-LRR")] = "NBS-LRR"
tx$cat2[tx$cat2 %in% c("NCR", "CRP0000-1030", "CRP1600-6250")] = "CRP"

fam = "CRP"
  tm = tx[tx$cat2 == fam,]
  tm1 = ddply(tm, .(clu), summarise, csize = length(clu))
  tm2 = ddply(tm1, .(csize), summarise, nc = length(csize), ng = sum(csize))
  tm3 = cbind(tm2, itv = sapply(tm2$csize, toString), stringsAsFactors = F)
  tm3$itv[tm3$csize>=11 & tm3$csize<=15] = '11-15'
  tm3$itv[tm3$csize>=16] = '16+'
  
sprintf("%s: %d, %.03f in tandem clusters\n", fam, nrow(tm), sum(tm2$ng[tm2$csize>1])/nrow(tm))
