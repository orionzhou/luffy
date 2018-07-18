require(plyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(RColorBrewer)
source('comp.fun.R')

dirw = file.path(Sys.getenv("misc3"), "comp.ortho")
diro = file.path(Sys.getenv("misc3"), "comp.stat")

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

ti = tst[tca$cat2 != 'TE',]

norgs = apply(ti[,orgs], 1, function(x) sum(!x %in% c('','-')))
table(norgs)

ti = cbind(ti, norg = norgs)

##### plot pan-proteome AFS
tab1 = table(ti$norg)
tab1 = tab1[names(tab1) != 1]
dt1 = data.frame(norg = as.numeric(names(tab1)), cnt = as.numeric(tab1), org = 'mixed', stringsAsFactors = F)

tis = ti[ti$norg == 1,]
oorg = apply(tis, 1, function(z) orgs[which(z[2:ncol(tis)] %in% c('syn', 'rbh'))])
tab2 = table(oorg)
dt2 = data.frame(norg = 1, cnt = as.numeric(tab2), org = names(tab2), stringsAsFactors = F)

to = rbind(dt1, dt2)

cols = c(brewer.pal(12, 'Set3'), brewer.pal(4, 'Set1'), 'gray30')
labs = orgs

to$org = factor(to$org, levels = c(orgs, 'mixed'))
to$norg = factor(to$norg, levels = sort(as.numeric(unique(to$norg))))
p1 = ggplot(to, aes(x = norg, y = cnt, fill = org, order = plyr:::desc(org))) +
  geom_bar(stat = 'identity', position = "stack", geom_params=list(width = 0.5)) +
  scale_fill_manual(name = "Accession-Specific", breaks = labs, labels = labs, values = cols, guide = guide_legend(ncol = 1, byrow = F, label.position = "right", direction = "vertical", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  scale_x_discrete(name = '# Sharing Accession') +
  scale_y_continuous(name = '# Gene Models', expand = c(0, 0), limits = c(0, 34000)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(legend.position = c(0.7, 0.7), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.6, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_text(size = 8, angle = 0), legend.text = element_text(size = 8, angle = 0)) +
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 90, hjust = 1))

#fp = sprintf("%s/71.pan.proteome.afs.pdf", diro)
#ggsave(p, filename = fp, width = 5, height = 5)


##### pan-proteome size
## run simulation
reps = 1:8
n_orgs = 1:16
tp = data.frame(rep = rep(reps, each = length(n_orgs)), 
  n_org = rep(n_orgs, length(rep)), core = NA, pan = NA)

for (rep in reps) {
  set.seed(rep * 100)
  orgs1 = c(tname, sample(qnames))
  tis = ti[,orgs1]
  tp$pan[tp$rep == rep & tp$n_org == 1] = sum(!tis[,tname] %in% c('', '-'))
  tp$core[tp$rep == rep & tp$n_org == 1] = sum(!tis[,tname] %in% c('', '-'))
  for (i in 2:length(orgs1)) {
    orgss = orgs1[1:i]
    norgs = apply(tis[,orgss], 1, function(z) sum(!z %in% c('', '-')))
    pan = sum(norgs >= 1)
    core = sum(norgs == length(orgss))
    tp$pan[tp$rep == rep & tp$n_org == i] = pan
    tp$core[tp$rep == rep & tp$n_org == i] = core
    cat(rep, orgs1[i], core, pan, "\n")
  }
}
fo = file.path(diro, "73.pan.proteome.size.tbl")
write.table(tp, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

## plot
fi = file.path(diro, "73.pan.proteome.size.tbl")
tp = read.table(fi, header = T, sep = "\t", as.is = T)

tp$rep = factor(tp$rep, levels = 1:max(tp$rep))
p2 = ggplot(tp) +
  geom_point(aes(x = n_org, y = pan), shape = 1, size = 1) +
  geom_point(aes(x = n_org, y = core), shape = 4, size = 1) +
  stat_smooth(aes(x = n_org, y = pan, col = 'a'), size = 0.3, se = F) +
  stat_smooth(aes(x = n_org, y = core, col = 'b'), size = 0.3, se = F) +
  scale_color_manual(name = "", labels = c('Pan-proteome', 'Core-proteome'), values = c("dodgerblue", "firebrick1")) +
  scale_x_continuous(name = '# Genomes Sequenced') +
  scale_y_continuous(name = '# Genes', expand = c(0, 0), limits = c(0, 110000)) + 
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.2, 'lines')) +
  theme(legend.position = c(0.15, 0.9), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA), legend.key.size = unit(1, 'lines'), legend.margin = unit(0, "line"), legend.title = element_blank(), legend.text = element_text(size = 8, angle = 0)) +
  theme(plot.margin = unit(c(0,1,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, color = "blue")) +
  theme(axis.text.y = element_text(size = 8, color = "grey", angle = 90, hjust  = 0.5))
  
#fp = sprintf("%s/73.pan.proteome.size.pdf", diro)
#ggsave(p, filename = fp, width = 5, height = 4)

fo = sprintf("%s/73.pan.proteome.pdf", diro)
pdf(file = fo, width = 10, height = 5, bg = 'transparent')
numrow = 1; numcol = 2
grid.newpage()
pushViewport(viewport(layout = grid.layout(numrow, numcol)))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

dco = data.frame(x = rep(1:numrow, each = numcol), y = rep(1:numcol, numrow), lab = LETTERS[1:(numrow*numcol)])
for (i in 1:nrow(dco)) {
  x = dco$x[i]; y = dco$y[i]; lab = dco$lab[i]
  grid.text(lab, x = 0, y = unit(1, 'npc'), just = c('left', 'top'), gp = gpar(col = "black", fontface = 2, fontsize = 20),
    vp = viewport(layout.pos.row = x, layout.pos.col = y))
}
dev.off()
