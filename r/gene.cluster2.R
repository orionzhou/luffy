require(grid)
require(plyr)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
source('Location.R')
source("comp.fun.R")

dirw = file.path(Sys.getenv("misc2"), "gene.cluster")
dir.create(dirw)

for (org in qnames_alpaca_comp) {
#org = "HM340.PJ"
f_gen = file.path(Sys.getenv("genome"), org, "51.gtb")
tg = read.table(f_gen, header = T, sep = "\t", as.is = T)[,c(1,3:5,16:17)]

### run Usearch and process output
f_fas = file.path(Sys.getenv("genome"), org, "51.fas")
cmd = sprintf("pro.cluster.py %s %s/61.usearch/%s 0.6", f_fas, dirw, org)
system(cmd)

fi = sprintf("%s/61.usearch/%s/32.tbl", dirw, org)
ti = read.table(fi, sep = "\t", header = T, as.is = T)
tc = merge(tg, ti, by = 'id', all.x = T)
y = which(is.na(tc$grp))
tc$grp[y] = seq(from = max(tc$grp, na.rm = T)+1, by = 1, length.out = length(y))

tt = tc[tc$cat2 != 'TE',]
tt = tt[order(tt$chr, tt$beg, tt$end),]
tt = cbind(idx = 1:nrow(tt), tt)
to = merge(tc[,c('id','grp')], tt[,c('idx','id')], by = 'id', all.x = T)

fo = sprintf("%s/61.usearch/%s.tbl", dirw, org)
write.table(to, fo, col.names = T, row.names = F, sep = "\t", quote = F)

### identify tandem clusters
fi = sprintf("%s/61.usearch/%s.tbl", dirw, org)
ti = read.table(fi, sep = "\t", header = T, as.is = T)
ti = ti[!is.na(ti$idx),]
ti = ti[order(ti$idx),]

aids = c(); bids = c()
wsize = 2
for (i in 1:nrow(ti)) {
  for (j in (i+1):(i+wsize)) {
    if(j > nrow(ti)) next
    if(ti$grp[i] == ti$grp[j]) {
      aids = c(aids, ti$id[i])
      bids = c(bids, ti$id[j])
    }
  }
}

tt = data.frame(aid = aids, bid = bids, stringsAsFactors = F)
write.table(tt, 'tmp.tbl', col.names = T, row.names = F, sep = "\t", quote = F)
system("graph.connComp.py tmp.tbl tmpo.tbl")
tc = read.table("tmpo.tbl", sep = "\t", header = T, as.is = T)
system(sprintf("rm %s %s", "tmp.tbl", "tmpo.tbl"))

to = merge(tg, tc, by = 'id', all.x = T)
fo = sprintf("%s/63.tandem/%s.tbl", dirw, org)
write.table(to, fo, col.names = T, row.names = F, sep = "\t", quote = F, na = '')

}

##### tandem array stats in  different genomes
to = data.frame()
for (org in qnames_alpaca_comp) {
fi = sprintf("%s/63.tandem/%s.tbl", dirw, org)
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti = ti[ti$cat2 != 'TE' & !is.na(ti$clu),]

gb = group_by(ti, clu)
tc = summarise(gb, csize = n())
ti = merge(ti, tc, by = 'clu')

brks = c(seq(1.5, 6.5, by = 1), Inf)
labs = c(2:7, '8+')
labs = factor(labs, levels = labs)

tx = ti
tx$cat2[tx$cat2 %in% c("CC-NBS-LRR", "TIR-NBS-LRR")] = "NBS-LRR"
fams = c("NBS-LRR", "F-box", "LRR-RLK", "NCR", "Unknown", "CRP0000-1030", "CRP1600-6250")
tx$cat2[! tx$cat2 %in% fams] = 'Pfam-Miscellaneous'

do = data.frame()
for (fam in c('NCR', 'NBS-LRR', 'LRR-RLK', "F-box", 'Pfam-Miscellaneous')) {
  tm = tx[tx$cat2 == fam,]
  tm2 = ddply(tm, .(csize), summarise, nc = length(csize), ng = sum(csize))
  tm3 = cbind(tm2, itv = sapply(tm2$csize, toString), stringsAsFactors = F)    
  tm3$itv[tm3$csize > brks[length(brks)-1]+1] = as.character(labs[length(labs)])
  tm4 = ddply(tm3, .(itv), summarise, nc = sum(nc), ng = sum(ng))
  tm5 = merge(data.frame(itv = labs), tm4, all.x = T)
  tm5$nc[is.na(tm5$nc)] = 0
  tm5$ng[is.na(tm5$ng)] = 0
  ds = data.frame(org = org, fam = fam, tm5, prop = tm5$ng/sum(tm5$ng), stringsAsFactors = T)
  do = rbind(do, ds)
}
to = rbind(to, do)
}
to$org = factor(to$org, levels = qnames_alpaca_comp)
to$itv = factor(to$itv, levels = labs)
x = to[to$itv == '2' & to$org == 'HM101',]
fams = as.character(x$fam[order(x$prop, decreasing = T)])
to$fam = factor(to$fam, levels = fams)

tw = ddply(to, .(org, itv), summarise, ng = sum(ng))
tw = reshape(tw, direction = 'wide', timevar = c('itv'), idvar = c('org'))
fo = file.path(dirw, "69.tbl")
#write.table(tw, fo, col.names = T, row.names = F, sep = "\t", quote = F)

cols = c(brewer.pal(3, "Blues"), brewer.pal(3, "Greens"), brewer.pal(3, "Oranges"), "Purple")

p1 = ggplot(to) + 
  geom_bar(aes(x = itv, y = ng, fill = org), stat = 'identity', position = 'dodge', width = 0.8) + 
  scale_x_discrete(name = 'Tandem Array Size') +
  scale_y_continuous(name = '# Tandem Arrays') +
  scale_fill_brewer(palette = "Paired") +
  facet_grid(fam ~ ., scales = 'free') +  
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(legend.position = "top", legend.title = element_blank(), legend.key.size = unit(0.6, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, angle = 0, colour = "blue", hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 0.5))

fo = sprintf("%s/69.pdf", dirw)
ggsave(p1, filename = fo, width = 6, height = 8)

to2 = ddply(to, .(org, itv), summarise, nc = sum(nc), ng = sum(ng))
to2$org = factor(to2$org, levels = qnames_alpaca_comp)
p1 = ggplot(to2) + 
  geom_bar(aes(x = itv, y = nc, fill = org), stat = 'identity', position = 'dodge', width = 0.8) + 
  scale_x_discrete(name = 'Tandem Array Size') +
  scale_y_continuous(name = '# Tandem Arrays') +
  scale_fill_manual(values = cols, guide = guide_legend(nrow = 3)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(legend.position = c(0.7, 0.9), legend.title = element_blank(), legend.key.size = unit(0.6, 'lines'), legend.background = element_rect(fill = 'white', size=0)) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, angle = 0, colour = "blue", hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 0.5))

fo = sprintf("%s/69.all.pdf", dirw)
ggsave(p1, filename = fo, width = 7, height = 5)


do = data.frame()
for (org in qnames_alpaca_comp) {
f_gen = file.path(Sys.getenv("genome"), org, "51.gtb")
tg = read.table(f_gen, header = T, sep = "\t", as.is = T)[,c(1,3:6,16:17)]

fi = sprintf("%s/63.tandem/%s.tbl", dirw, org)
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti = ti[ti$cat2 != 'TE' & !is.na(ti$clu),]

to = merge(tg, ti[,c('id','clu')], by = 'id')
to = cbind(org = org, to)
to$clu = sprintf("%s_clus_%d", org, to$clu)
do = rbind(do, to)

pct = nrow(ti)/sum(tg$cat2 != "TE")
cat(sprintf("%10s: %5d genes, %d clusters | %d (%.03f) in tandem arrays\n", org, nrow(tg), length(unique(ti$clu)), nrow(ti), pct))
}
fo = file.path(dirw, "69.2.tbl")
write.table(do, fo, col.names = T, row.names = F, sep = "\t", quote = F)
