require(plyr)
require(dplyr)
require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(xlsx)
require(seqinr)
require(RColorBrewer)
require(ape)
source("Location.R")
source("comp.fun.R")

diro = file.path(Sys.getenv('misc3'), 'alpaca')
setwd(diro)

qnames = qnames_alpaca_comp

tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

## create a config file
fps = sprintf("%s/%s/51.fas", Sys.getenv('genome'), qnames)
to = data.frame(org = qnames, fp = fps, stringsAsFactors = F)
fo = file.path(diro, "31.conf.csv")
#write.table(to, fo, sep = ",", row.names = F, col.names = F, quote = F)
#merge.fas.py 31.conf.csv 32.pro.fas

## read all genes
tg = data.frame()
for (qname in qnames) {
  dirg = sprintf("%s/%s", Sys.getenv("genome"), qname)
  fg = file.path(dirg, "51.gtb")
  tg1 = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1,3:6,16:17)]
  tg = rbind(tg, cbind(org = qname, tg1))
}
colnames(tg)[c(2,7:8)] = c('gid','fam','sfam')

## write gene ids 
sfams = c('CRP0355', 'CRP3710', 'CRP4180')
opts = c('alps', 'alpc')

for (sfam in sfams) {
  tgs = tg[tg$sfam == sfam,]
  tgs = cbind(tgs, id = sprintf("%s|%s", tgs$org, tgs$gid))
  ids1 = tgs$id[tgs$org %in% c("HM101", "HM056", "HM034", "HM340")]
  ids2 = tgs$id[tgs$org %in% c("HM101", "HM056.AC", "HM034.AC", "HM340.AC")]
  fo1 = sprintf("%s/34.sfams/%s.alps.bed", diro, sfam)
  fo2 = sprintf("%s/34.sfams/%s.alpc.bed", diro, sfam)
  write.table(ids1, fo1, sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(ids2, fo2, sep = "\t", row.names = F, col.names = F, quote = F)
}

### run alpaca.genefam.py
qnames1 = c("HM101","HM056","HM034","HM340")
qnames2 = c("HM101","HM056.AC","HM034.AC","HM340.AC")
cols_raw = brewer.pal(10, "Paired")
cols_map1 = cols_raw[c(10,1,3,5)]
cols_map2 = cols_raw[c(10,2,4,6)]
names(cols_map1) = qnames1
names(cols_map2) = qnames2


sfam = sfams[1]

  ft1 = sprintf("%s/34.sfams/%s.alps.phy.nwk", diro, sfam)
  tree1 = read.tree(ft1)
  res = strsplit(tree1$tip.label, '[|]')
  oids = sapply(res, "[", 1)
  chrs = sapply(strsplit(sapply(res, "[", 2), "_"), "[", 2)
  tree1$tip.label = chrs
  font1 = rep(1, length(oids))
  tip.cols1 = cols_map1[oids]

  ft2 = sprintf("%s/34.sfams/%s.alpc.phy.nwk", diro, sfam)
  tree2 = read.tree(ft2)
  res = strsplit(tree2$tip.label, '[|]')
  oids = sapply(res, "[", 1)
  chrs = sapply(strsplit(sapply(res, "[", 2), "_"), "[", 2)
  tree2$tip.label = chrs
  font2 = rep(1, length(oids))
  tip.cols2 = cols_map2[oids]
  
  fo = sprintf("%s/34.sfams/%s.pdf", diro, sfam)
  pdf(fo, width = 12, height = 10)
  par(mfrow = c(1,2))
  
  plot(tree1, type = 'phylogram', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols1, label.offset = 0.11, font = font1,
    no.margin = T, cex = 0.6)
  tiplabels(pch = 22, frame = 'none', adj = 0.55, bg = tip.cols1)
  add.scale.bar(x = 0.02, y = tree1$Nnode*0.9 , lcol = 'black')
  legend(0.02, tree1$Nnode*0.8, legend=qnames1, fill=cols_map1, bty='n', cex=0.8)
  text(2, tree1$Nnode*1.01, "ALLPATHS", cex=1.5)
  
  plot(tree2, type = 'phylogram', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols2, label.offset = 0.11, font = font2,
    no.margin = T, cex = 0.6)
  tiplabels(pch = 22, frame = 'none', adj = 0.55, bg = tip.cols2)
  add.scale.bar(x = 0.02, y = tree2$Nnode*0.9 , lcol = 'black')
  legend(0.02, tree2$Nnode*0.8, legend=qnames2, fill=cols_map2, bty='n', cex=0.8)
  text(2, tree2$Nnode*1.01, "ALPACA", cex=1.5)

  dev.off()


sfam = sfams[2]

  ft1 = sprintf("%s/34.sfams/%s.alps.phy.nwk", diro, sfam)
  tree1 = read.tree(ft1)
  res = strsplit(tree1$tip.label, '[|]')
  oids = sapply(res, "[", 1)
  chrs = sapply(strsplit(sapply(res, "[", 2), "_"), "[", 2)
  tree1$tip.label = chrs
  font1 = rep(1, length(oids))
  tip.cols1 = cols_map1[oids]

  ft2 = sprintf("%s/34.sfams/%s.alpc.phy.nwk", diro, sfam)
  tree2 = read.tree(ft2)
  res = strsplit(tree2$tip.label, '[|]')
  oids = sapply(res, "[", 1)
  chrs = sapply(strsplit(sapply(res, "[", 2), "_"), "[", 2)
  tree2$tip.label = chrs
  font2 = rep(1, length(oids))
  tip.cols2 = cols_map2[oids]
  
  fo = sprintf("%s/34.sfams/%s.pdf", diro, sfam)
  pdf(fo, width = 12, height = 10)
  par(mfrow = c(1,2))
  
  plot(tree1, type = 'phylogram', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols1, label.offset = 0.11, font = font1,
    no.margin = T, cex = 0.6)
  tiplabels(pch = 22, frame = 'none', adj = 0.55, bg = tip.cols1)
  add.scale.bar(x = 0.02, y = tree1$Nnode*0.95 , lcol = 'black')
  legend(0.02, tree1$Nnode*0.9, legend=qnames1, fill=cols_map1, bty='n', cex=0.8)
  text(1, tree1$Nnode*1.02, "ALLPATHS", cex=1.5)
  
  plot(tree2, type = 'phylogram', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols2, label.offset = 0.11, font = font2,
    no.margin = T, cex = 0.6)
  tiplabels(pch = 22, frame = 'none', adj = 0.55, bg = tip.cols2)
  add.scale.bar(x = 0.02, y = tree2$Nnode*0.9 , lcol = 'black')
  legend(0.02, tree2$Nnode*0.8, legend=qnames2, fill=cols_map2, bty='n', cex=0.8)
  text(1, tree2$Nnode*1.01, "ALPACA", cex=1.5)

  dev.off()


sfam = sfams[3]
  ft1 = sprintf("%s/34.sfams/%s.alps.phy.nwk", diro, sfam)
  tree1 = read.tree(ft1)
  res = strsplit(tree1$tip.label, '[|]')
  oids = sapply(res, "[", 1)
  chrs = sapply(strsplit(sapply(res, "[", 2), "_"), "[", 2)
  tree1$tip.label = chrs
  font1 = rep(1, length(oids))
  tip.cols1 = cols_map1[oids]

  ft2 = sprintf("%s/34.sfams/%s.alpc.phy.nwk", diro, sfam)
  tree2 = read.tree(ft2)
  res = strsplit(tree2$tip.label, '[|]')
  oids = sapply(res, "[", 1)
  chrs = sapply(strsplit(sapply(res, "[", 2), "_"), "[", 2)
  tree2$tip.label = chrs
  font2 = rep(1, length(oids))
  tip.cols2 = cols_map2[oids]
  
  fo = sprintf("%s/34.sfams/%s.pdf", diro, sfam)
  pdf(fo, width = 12, height = 10)
  par(mfrow = c(1,2))
  
  plot(tree1, type = 'phylogram', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols1, label.offset = 0.06, font = font1,
    no.margin = T, cex = 0.6)
  tiplabels(pch = 22, frame = 'none', adj = 0.52, bg = tip.cols1)
  add.scale.bar(x = 0.02, y = tree1$Nnode*0.9 , lcol = 'black')
  legend(0.02, tree1$Nnode*0.8, legend=qnames1, fill=cols_map1, bty='n', cex=0.8)
  text(2, tree1$Nnode*1.01, "ALLPATHS", cex=1.5)
  
  plot(tree2, type = 'phylogram', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols2, label.offset = 0.06, font = font2,
    no.margin = T, cex = 0.6)
  tiplabels(pch = 22, frame = 'none', adj = 0.52, bg = tip.cols2)
  add.scale.bar(x = 0.02, y = tree2$Nnode*0.9 , lcol = 'black')
  legend(0.02, tree2$Nnode*0.8, legend=qnames2, fill=cols_map2, bty='n', cex=0.8)
  text(2, tree2$Nnode*1.01, "ALPACA", cex=1.5)
  
  dev.off()
  
  png(filename = fo1, width = 600, height = 800, units = 'px')
  plot(tree, type = 'phylogram', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols, label.offset = 0.06, font = font,
    no.margin = T, cex = 1)
  tiplabels(pch = 22, frame = 'none', adj = 0.52, bg = tip.cols)
 # nodelabels(pch = 22, bg = node.bg)
  add.scale.bar(x = 0.02, y = tree$Nnode*0.9 , lcol = 'black')
