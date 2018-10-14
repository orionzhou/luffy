require(grid)
require(plyr)
require(dplyr)
require(ggplot2)
source('Location.R')
source("comp.fun.R")

dirw = file.path(Sys.getenv("misc2"), "gene.cluster")
dir.create(dirw)

org = "HM101"

f_gen = file.path(Sys.getenv("genome"), org, "51.gtb")
tg = read.table(f_gen, header = T, sep = "\t", as.is = T)[,c(1,3:5,16:17)]

### process mcl output
fi = sprintf("%s/04.mcl/%s.mcl", dirw, org)
ti = read.table(fi, sep = "-", header = F, as.is = T)
ti = cbind(grp = 1:nrow(ti), ti)
x = apply(ti, 1, kkk <- function(rw) { 
  data.frame(grp = as.numeric(rw[1]), id = unlist(strsplit(rw[2], "\t")), stringsAsFactors = F)
})
require(data.table)
tc = rbindlist(x)
tc = merge(tg, tc, by = 'id', all.x = T)
y = which(is.na(tc$grp))
tc$grp[y] = seq(from = max(tc$grp, na.rm = T)+1, by = 1, length.out = length(y))

tt = tc[tc$cat2 != 'TE',]
tt = tt[order(tt$chr, tt$beg, tt$end),]
tt = cbind(idx = 1:nrow(tt), tt)
to = merge(tc[,c('id','grp')], tt[,c('idx','id')], by = 'id', all.x = T)

fo = sprintf("%s/08.grp/%s.tbl", dirw, org)
write.table(to, fo, col.names = T, row.names = F, sep = "\t", quote = F)


### identify tandem clusters
fi = sprintf("%s/08.grp/%s.tbl", dirw, org)
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
fo = sprintf("%s/11.tandem/%s.tbl", dirw, org)
write.table(to, fo, col.names = T, row.names = F, sep = "\t", quote = F, na = '')

### cluster distribution for different families
fi = sprintf("%s/11.tandem/%s.tbl", dirw, org)
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti = ti[ti$cat2 != 'TE',]
idx_single = which(is.na(ti$clu))
ti$clu[idx_single] = seq(max(ti$clu, na.rm = T)+1, by = 1, length.out = length(idx_single))

gb = group_by(ti, clu)
tc = summarise(gb, csize = n())
ti = merge(ti, tc, by = 'clu')

brks = c(seq(0.5, 10.5, by = 1), 15.5, Inf)
labs = c(1:10, '11-15', '16+')
labs = factor(labs, levels = labs)

to = ti
to$cat2[to$cat2 %in% c("CC-NBS-LRR", "TIR-NBS-LRR")] = "NBS-LRR"
fams = c("NBS-LRR", "F-box", "LRR-RLK", "NCR", "Unknown", "CRP0000-1030", "CRP1600-6250")
to$cat2[! to$cat2 %in% fams] = 'Pfam-Miscellaneous'

do = data.frame()
for (fam in c('NCR', 'NBS-LRR', 'LRR-RLK', "F-box", 'Pfam-Miscellaneous')) {
  tm = to[to$cat2 == fam,]
  tm2 = ddply(tm, .(csize), summarise, nc = length(csize), ng = sum(csize))
  tm3 = cbind(tm2, itv = sapply(tm2$csize, toString), stringsAsFactors = F)
  tm3$itv[tm3$csize>=11 & tm3$csize<=15] = '11-15'
  tm3$itv[tm3$csize>=16] = '16+'
  tm4 = ddply(tm3, .(itv), summarise, nc = sum(nc), ng = sum(ng))
  tm5 = merge(data.frame(itv = labs), tm4, all.x = T)
  tm5$nc[is.na(tm5$nc)] = 0; tm5$ng[is.na(tm5$ng)] = 0
  ds = data.frame(org = org, fam = fam, tm5, prop = tm5$ng/sum(tm5$ng), stringsAsFactors = T)
  do = rbind(do, ds)
}

x = do[do$itv == '1',]
fams = as.character(x$fam[order(x$prop, decreasing = T)])
do$fam = factor(do$fam, levels = fams)
do$itv = factor(do$itv, levels = labs)

p1 = ggplot(do[do$fam == 'Pfam-Miscellaneous',]) + 
  geom_bar(aes(x = itv, y = ng), stat = 'identity', geom_params = list(width = 0.8)) + 
  scale_x_discrete(name = 'Tandem array size') +
  scale_y_continuous(name = '# Genes') +
  facet_wrap(~ fam, scales = 'fixed', nrow = 1) +  
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
#  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, angle = 60, colour = "blue", hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 0.5))
p2 = ggplot(do[do$fam != 'Pfam-Miscellaneous',]) + 
  geom_bar(aes(x = itv, y = ng), stat = 'identity', geom_params = list(width = 0.8)) + 
  scale_x_discrete(name = 'Tandem array size') +
  scale_y_continuous(name = '# Genes') +
  facet_wrap(~ fam, scales = 'fixed', nrow = 1) +  
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(plot.margin = unit(c(0,0.5,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, angle = 60, colour = "blue", hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 0.5))

fo = sprintf("%s/12.tandem.plot/%s.pdf", dirw, org)
pdf(file = fo, width = 8, height = 5, bg = 'transparent')
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3, width = c(1,0.95,1), heights = c(1,1))))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:3)) 
dev.off()


##### tandem array stats in  different genomes
to = data.frame()
orgs = c(tname, "HM056", "HM056.AC", "HM034", "HM034.AC", "HM340", "HM340.AC")
for (org in orgs) {
fi = sprintf("%s/11.tandem/%s.tbl", dirw, org)
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti = ti[ti$cat2 != 'TE',]
idx_single = which(is.na(ti$clu))
ti$clu[idx_single] = seq(max(ti$clu, na.rm = T)+1, by = 1, length.out = length(idx_single))

gb = group_by(ti, clu)
tc = summarise(gb, csize = n())
ti = merge(ti, tc, by = 'clu')

brks = c(seq(0.5, 10.5, by = 1), 15.5, Inf)
labs = c(1:10, '11-15', '16+')
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
  tm3$itv[tm3$csize>=11 & tm3$csize<=15] = '11-15'
  tm3$itv[tm3$csize>=16] = '16+'
  tm4 = ddply(tm3, .(itv), summarise, nc = sum(nc), ng = sum(ng))
  ds = data.frame(org = org, fam = fam, tm4, prop = tm4$ng/sum(tm4$ng), stringsAsFactors = T)
  do = rbind(do, ds)
}
to = rbind(to, do)
}
to$org = factor(to$org, levels = orgs)
to$itv = factor(to$itv, levels = labs)
x = to[to$itv == '1' & to$org == 'HM101',]
fams = as.character(x$fam[order(x$prop, decreasing = T)])
to$fam = factor(to$fam, levels = fams)

fo1 = file.path(dirw, "50.tbl")
write.table(to, fo1, col.names = T, row.names = F, sep = "\t", quote = F)

tw = ddply(to, .(org, itv), summarise, ng = sum(ng))
tw = reshape(tw, direction = 'wide', timevar = c('itv'), idvar = c('org'))
fo = file.path(dirw, "51.tbl")
write.table(tw, fo, col.names = T, row.names = F, sep = "\t", quote = F)

p1 = ggplot(to) + 
  geom_bar(aes(x = itv, y = prop), stat = 'identity', geom_params = list(width = 0.7)) + 
  scale_x_discrete(name = 'Tandem array size') +
  scale_y_continuous(name = 'Proportion in family') +
  facet_grid(org ~ fam)+#, scales = 'free', nrow = 2) +  
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, angle = 60, colour = "blue", hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 0.5))

fo = sprintf("%s/51.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 6)

for (org in orgs) {
f_gen = file.path(Sys.getenv("genome"), org, "51.gtb")
tg = read.table(f_gen, header = T, sep = "\t", as.is = T)[,c(1,3:5,15:17)]

fi = sprintf("%s/11_tandem/%s.tbl", dirw, org)
ti = read.table(fi, header = F, sep = "\t", as.is = T)
colnames(ti) = c('clu', 'id')

pct = nrow(ti)/nrow(tg)
cat(sprintf("%10s: %5d genes, %d clusters | %d (%.03f) in tandem arrays\n", org, nrow(tg), length(unique(ti$clu)), nrow(ti), pct))
}
