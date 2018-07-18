require(plyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(reshape2)
require(RColorBrewer)
require(ape)
require(bios2mds)
require(gridBase)
require(colorRamps)
source("comp.fun.R")


dirw = file.path(Sys.getenv("misc3"), "comp.ortho.hm")
diro = file.path(Sys.getenv("misc3"), "comp.genefam")

fg = file.path(tcfg$gdir, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1,16:17)]

fi = file.path(dirw, "29.mpd.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)

sum(tg$id != ti$id)
tm = cbind(tg[,2:3], mpd = ti$mpd)
colnames(tm)[1:2] = c('fam', 'sfam')
tm = tm[!is.na(tm$mpd),]

fo = file.path(diro, "12.mpd.tbl")
write.table(tm, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')


tis = ddply(tm, .(fam), summarise, cnt = length(fam), q25 = quantile(mpd, 0.25), q50 = quantile(mpd, 0.5), q75 = quantile(mpd, 0.75))
tis = tis[order(tis$q50, decreasing = T),]
fams = tis$fam

tis$fam = factor(tis$fam, levels = fams)
p1 = ggplot(tis) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', position = 'dodge', width = 0.6) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01), labels = sprintf("%s | %5d", tis$fam, tis$cnt)) +
  scale_y_continuous(name = 'Mean Pariwise Protein Distance', expand = c(0,0), limits = c(0, 0.3)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = file.path(diro, "12.mpd.pdf")
ggsave(p1, filename = fp, width = 5, height = 8)



### obsolete
dirw = file.path(Sys.getenv("misc3"), "comp.ortho")
diro = file.path(Sys.getenv("misc3"), "comp.genefam")

#####
fi = file.path(dirw, "21.gid.tbl")
tid = read.table(fi, header = T, sep = "\t", as.is = T)
fi = file.path(dirw, "21.loc.tbl")
tlo = read.table(fi, header = T, sep = "\t", as.is = T)
fi = file.path(dirw, "21.sta.tbl")
tst = read.table(fi, header = T, sep = "\t", as.is = T)
fi = file.path(dirw, "22.cat.tbl")
tca = read.table(fi, header = T, sep = "\t", as.is = T)
fi = file.path(dirw, "28.dist.tbl")
tds = read.table(fi, header = T, sep = "\t", as.is = T)

norgs = apply(tst[,qnames_15], 1, function(x) sum(!x %in% c('','-')))
table(norgs)

##### gene fam mpd distribution 
orgs_in = c(tname, qnames_12)
  comps = c()
  ids = orgs_in
  for (i in 1:(length(ids)-1)) {
    id1 = ids[i]
    for (j in (i+1):length(ids)) {
      id2 = ids[j]
      comps = c(comps, sprintf("%s.%s", id1, id2))
    }
  }

#### cross-bar plot showing quantile
norgs = apply(tst[,qnames_15], 1, function(x) sum(!x %in% c('','-')))
idxs = which(norgs > 10)
mpd = apply(tds[idxs,comps], 1, mean, na.rm = T)
idxs = idxs[!is.na(mpd)]
mpd = mpd[!is.na(mpd)]
tu = data.frame(fam = tca$cat2[idxs], sfam = tca$cat3[idxs], mpd = mpd, stringsAsFactors = F)

fo = file.path(diro, "12.mpd.tbl")
write.table(tu, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

tis = ddply(tu, .(fam), summarise, cnt = length(fam), q25 = quantile(mpd, 0.25), q50 = quantile(mpd, 0.5), q75 = quantile(mpd, 0.75))
tis = tis[order(tis$q50, decreasing = T),]
fams = tis$fam

tis$fam = factor(tis$fam, levels = fams)
p1 = ggplot(tis) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', position = 'dodge', width = 0.6) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01), labels = sprintf("%s | %5d", tis$fam, tis$cnt)) +
  scale_y_continuous(name = 'Mean Pariwise Protein Distance', expand = c(0,0), limits = c(0, 0.2)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = file.path(diro, "12.mpd.pdf")
ggsave(p1, filename = fp, width = 5, height = 8)

