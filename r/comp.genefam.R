require(plyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(RColorBrewer)
require(gtable)
source('comp.fun.R')

dirw = file.path(Sys.getenv("misc3"), "comp.genefam")

ffam = "/home/youngn/zhoux379/data/db/pfam/genefam.tsv"
tfam = read.table(ffam, header = T, sep = "\t", as.is = T)[,1:2]

tg = read.table(tcfg$gene, sep = "\t", header = F, as.is = T)
colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
tg = tg[tg$type == 'cds',]

##### genefam three-panel plot: thetaPi + largeEff + mpd + cnv
# pi
fi = file.path(dirw, "03.pi.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
to = rename_genefam(ti)

tp = ddply(to, .(fam), summarise, cnt = length(fam), q25 = quantile(pi, 0.25, na.rm = T), q50 = median(pi, na.rm = T), q75 = quantile(pi, 0.75, na.rm = T))
famsf = as.character(tp$fam[order(tp$q50, decreasing = T)])
tp$fam = factor(tp$fam, levels = famsf)
labs = sprintf("%s (%d)", famsf, tp$cnt[match(famsf, tp$fam)])

p1 = ggplot(tp) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', position = 'dodge', width = 0.6) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01), labels = labs) +
  scale_y_continuous(name = 'ThetaPi', expand = c(0, 0), limits = c(0, 0.02)) +
  theme_bw() +
  theme(plot.margin = unit(c(1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

# largeeff
fi = file.path(dirw, "07.largeeff.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
to = rename_genefam(ti)

tp = ddply(to, .(fam, eff), summarise, prop = sum(cnt > 0) / length(cnt))
tp$fam = factor(tp$fam, levels = famsf)

p2 = ggplot(tp) +
  geom_bar(aes(x = fam, y = prop, fill = eff),
    stat = 'identity', position = 'stack', width = 0.6) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01)) +
  scale_y_continuous(name = 'Proportion', expand = c(0, 0), limits = c(0, 1)) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  theme(legend.position = c(0.7, 0.8), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.7, 'lines'), legend.title = element_blank(), legend.text = element_text(size = 8, angle = 0)) +
  theme(plot.margin = unit(c(1,0.1,0.1,1.1), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
#  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))
  theme(axis.text.y = element_blank())

# mpd
fi = file.path(dirw, "12.mpd.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
to = rename_genefam(ti)

tp = ddply(to, .(fam), summarise, cnt = length(fam), q25 = quantile(mpd, 0.25), q50 = quantile(mpd, 0.5), q75 = quantile(mpd, 0.75))
tp = rbind(tp, data.frame(fam='LRR',cnt=0,q25=NA,q50=NA,q75=NA))
tp$fam = factor(tp$fam, levels = famsf)
labs = sprintf("%s (%d)", famsf, tp$cnt[match(famsf, tp$fam)])

p3 = ggplot(tp) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', position = 'dodge', width = 0.6) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01), labels = labs) +
  scale_y_continuous(name = 'Mean Pariwise Protein Distance', expand = c(0,0), limits = c(0, 0.23)) +
  theme_bw() +
  theme(plot.margin = unit(c(1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

# cnv
fi = file.path(dirw, "25.cnv.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
to = rename_genefam(ti)

tp = ddply(to, .(fam), summarise, cnt = length(fam), q25=quantile(cv,0.25,na.rm=T), q50=quantile(cv,0.5,na.rm=T), q75=quantile(cv,0.75,na.rm=T), avg=mean(cv,na.rm=T), std=sd(cv,na.rm=T))
tp$fam = factor(tp$fam, levels = famsf)
labs = sprintf("(%d)", tp$cnt[match(famsf, tp$fam)])

p4 = ggplot(tp) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', position = 'dodge', width = 0.6) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01), labels=labs) +
  scale_y_continuous(name = 'C.V. of Gene # in Ortholog Group', expand = c(0.01, 0.01)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0,'lines')) +
  theme(plot.margin = unit(c(1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

# multi-panel plot 
fo = file.path(dirw, "49.genefam.pdf")
numrow = 2; hts = c(3, 3)
numcol = 2; wds = c(4, 3)
pdf(file = fo, width = sum(wds), height = sum(hts), bg = 'transparent')
#tiff(file = fo, width = sum(wds), height = 4, units = 'in', bg = 'white')
grid.newpage()
pushViewport(viewport(layout = grid.layout(numrow, numcol, width = wds, heights = hts)))

ps = list(p1, p2, p3, p4)
xs = c(1,1,2,2)
ys = c(1,2,1,2)
ss = LETTERS[1:length(xs)]
for (i in 1:length(xs)) {
  x = xs[i]; y = ys[i]; lab = ss[i]
  print(ps[[i]], vp = viewport(layout.pos.row = x, layout.pos.col = y))
  grid.text(lab, x = 0, y = unit(1, 'npc'), just = c('left', 'top'), gp = gpar(col = "black", fontface = 2, fontsize = 20),
    vp = viewport(layout.pos.row = x, layout.pos.col = y))
}
dev.off()


# protein length (obsolete)
gb = dplyr::group_by(tg, id)
to = dplyr::summarise(gb, fam = fam[1], len = sum(end-beg+1))
to = rename_genefam(to, fams)

tp = ddply(to, .(fam), summarise, q25 = quantile(len, 0.25), q50 = quantile(len, 0.5), q75 = quantile(len, 0.75))
tp = tp[tp$fam %in% famsf,]
tp$fam = factor(tp$fam, levels = famsf)
#labs = sprintf("%s (%d)", famsf, tp$cnt[match(famsf, tp$fam)])

p5 = ggplot(tp) +
  geom_crossbar(aes(x = fam, y = q50/3, ymin = q25/3, ymax = q75/3),
    stat = 'identity', position = 'dodge', width = 0.6) +
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01)) +
  scale_y_continuous(name = 'Protein Length', expand = c(0, 0), limits = c(0, 1300)) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  theme(plot.margin = unit(c(1,1,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_blank())
 