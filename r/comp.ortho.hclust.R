require(ape)
require(plyr)
require(ggplot2)
require(reshape2)
require(Biostrings)
require(colorRamps)
require(gridBase)
source("Align.R")
source("Location.R")
source("comp.fun.R")

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
norgs = apply(tst[,orgs], 1, function(x) sum(! x%in% c('', '-')))

  ids = orgs
  xids = c(); yids = c()
  for (i in 1:(length(ids)-1)) {
    id1 = ids[i]
    for (j in (i+1):length(ids)) {
      id2 = ids[j]
      xids = c(xids, id1)
      yids = c(yids, id2)
    }
  }
lidx = list()
for (id in ids) {
  lidx[[id]] = which(xids == id | yids == id)
}

##### CRP ortho-map hclust
idxs = which(tca$cat2 %in% c("CRP0000-1030", "NCR", "CRP1600-6250") & norgs >= 1)
ti = matrix(NA, nrow = length(idxs), ncol = length(orgs), dimnames = list(NULL, orgs))

for (i in 1:length(idxs)) {
  idx = idxs[i]
  for (org in orgs) {
    sta = tst[idx, org]
    if(sta == '') {
      score = NA
    } else if(sta %in% c('-')) {
      score = 0
    } else {
      score = 1 - mean(as.numeric(tds[idx, lidx[[org]]]), na.rm = T)
      if(is.na(score)) {score = 1}
    }
    ti[i, org] = score
  }
}

### hclust
cordist <- function(x) as.dist(1-cor(t(x), method="pearson", use = "pairwise.complete.obs"))
r.dist <- cordist(t(ti))
r.hc <- hclust(r.dist, method='ward.D')
tree <- as.phylo(r.hc)

tree = rotate(tree, 19)
tree = rotate(tree, 20)
tree = reorder.phylo(tree, order = 'postorder')
plot(tree)

tw = cbind(i = 1:length(idxs), tlo[idxs, c('chr','pos')], fam = tca$cat3[idxs], fam2 = tca$cat2[idxs])
#tw$fam = factor(tw$fam, levels = c("CRP0000-1030", "NCR", "CRP1600-6250"))
tw = tw[order(tw$fam, tw$chr, tw$pos), ]
tf = cbind(x = 1:nrow(tw), tw)
tx = ddply(tf, .(fam2), summarise, beg = min(x), end = max(x), mid = as.integer(mean(x)))

#ff = '/home/youngn/zhoup/Data/misc2/genefam/crp.info'
#tf = read.table(ff, sep = "\t", as.is = T, header = T)
#tf = merge(tw[,c('cat3','x')], tf, by.x = 'cat3', by.y = 'id')

#tx = ddply(tw, .(chr), summarise, beg = min(x), end = max(x))
tx = tx[tx$end > tx$beg,]

tp = cbind(tf, ti[tf$i,])
tp = melt(tp[,c('x','fam',orgs)], id.vars = c('x', 'fam'), variable.name = 'org', value.name = 'score')

tp$org = factor(tp$org, levels = rev(tree$tip.label))

p = ggplot(tp, aes(x = x, y = org, fill = score)) +
  geom_tile(stat = 'identity', position = "identity") + 
  scale_fill_gradientn(colours = matlab.like2(20)) +
#  labs(fill = "Mean Pariwise AA distance") +
  scale_x_continuous(name = '', limits = c(0, max(tp$x)+1), expand=c(0, 0), breaks = floor((tx$beg+tx$end)/2), labels = as.character(tx$fam)) +
  scale_y_discrete(name = '') +
  theme(legend.position = "right", legend.key.size = unit(1, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.2, 'lines')) +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 8, color = "brown")) +
  theme(axis.text.y = element_text(size = 9, color = "black")) +
  annotate('segment', x = tx$beg, xend = tx$end, y = 0.3, yend = 0.3, size = 0.6, color = 'darkgreen') +
  annotate('segment', x = tx$beg, xend = tx$beg, y = 0.2, yend = 0.4, size = 0.5, color = 'darkgreen') +
  annotate('segment', x = tx$end, xend = tx$end, y = 0.2, yend = 0.4, size = 0.5, color = 'darkgreen')

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

fp = sprintf("%s/81.crp.pdf", diro)
pdf(file = fp, width = 10, height = 5, bg = 'transparent')
plot.new()

vl <- viewport(x = 0, y = 0, 
  width = unit(0.1, 'npc'), height = unit(1, 'npc'), 
  just = c('left', 'bottom'), name = 'left')
pushViewport(vl)
par(new = T, fig = gridFIG(), mar=c(1.3,0.4,0.2,0))
plot.phylo(tree, no.margin = F, show.tip.label = F)
upViewport()

vr <- viewport(x = unit(1, 'npc'), y = 0, 
  width = unit(0.9, 'npc'), height = unit(1, 'npc'), 
  just = c('right', 'bottom'), name = 'right')
pushViewport(vr)
grid.draw(gt)
upViewport()

dev.off()


##### plot NBS-LRR ortholog map
idxs = which(tca$cat2 %in% c('TIR-NBS-LRR', 'CC-NBS-LRR') & norgs >= 1)
#idxs = which(tca$cat2 %in% c('2OG-FeII_Oxy', 'Peroxidase', 'Cytochrome', 'Hydrolase', 'Peptidase'))

ti = matrix(NA, nrow = length(idxs), ncol = length(orgs), dimnames = list(NULL, orgs))

for (i in 1:length(idxs)) {
  idx = idxs[i]
  for (org in orgs) {
    sta = tst[idx, org]
    if(sta == '') {
      score = NA
    } else if(sta %in% c('-')) {
      score = 0
    } else {
      score = 1 - mean(as.numeric(tds[idx, lidx[[org]]]), na.rm = T)
      if(is.na(score)) {score = 1}
    }
    ti[i, org] = score
  }
}

### hclust
cordist <- function(x) as.dist(1-cor(t(x), method="pearson", use = "pairwise.complete.obs"))
r.dist <- cordist(t(ti))
r.hc <- hclust(r.dist, method='ward.D')
tree = as.phylo(r.hc)

tw = cbind(i = 1:length(idxs), tlo[idxs, c('chr','pos')], fam = tca$cat3[idxs], fam2 = tca$cat2[idxs])
#tw$fam = factor(tw$fam, levels = c("CRP0000-1030", "NCR", "CRP1600-6250"))
tw = tw[order(tw$fam, tw$chr, tw$pos), ]
tf = cbind(x = 1:nrow(tw), tw)
tx = ddply(tf, .(fam2), summarise, beg = min(x), end = max(x), mid = as.integer(mean(x)))
tx = tx[tx$end > tx$beg,]

tp = cbind(tf, ti[tf$i,])
tp = melt(tp[,c('x','fam',orgs)], id.vars = c('x', 'fam'), variable.name = 'org', value.name = 'score')

tp$org = factor(tp$org, levels = rev(tree$tip.label))

p = ggplot(tp, aes(x = x, y = org, fill = score)) +
  geom_tile(stat = 'identity', position = "identity") + 
  scale_fill_gradientn(colours = matlab.like2(20)) +
#  labs(fill = "Mean Pariwise AA distance") +
  scale_x_continuous(name = '', limits = c(0, max(tp$x)+1), expand=c(0, 0), breaks = floor((tx$beg+tx$end)/2), labels = tx$fam) +
  scale_y_discrete(name = '') +
  theme(legend.position = "right", legend.key.size = unit(1, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.2, 'lines')) +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 8, color = "brown")) +
  theme(axis.text.y = element_text(size = 9, color = "black")) +
  annotate('segment', x = tx$beg, xend = tx$end, y = 0.3, yend = 0.3, size = 0.6, color = 'darkgreen') +
  annotate('segment', x = tx$beg, xend = tx$beg, y = 0.2, yend = 0.4, size = 0.1, color = 'darkgreen') +
  annotate('segment', x = tx$end, xend = tx$end, y = 0.2, yend = 0.4, size = 0.1, color = 'darkgreen')

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

fp = sprintf("%s/81.nbs.pdf", diro)
#fp = sprintf("%s/81.misc.pdf", diro)
pdf(file = fp, width = 10, height = 5, bg = 'transparent')
plot.new()

vl <- viewport(x = 0, y = 0, 
  width = unit(0.1, 'npc'), height = unit(1, 'npc'), 
  just = c('left', 'bottom'), name = 'left')
pushViewport(vl)
par(new = T, fig = gridFIG(), mar=c(2,0.4,0.2,0))
plot.phylo(tree, no.margin = F, show.tip.label = F)
upViewport()

vr <- viewport(x = unit(1, 'npc'), y = 0, 
  width = unit(0.9, 'npc'), height = unit(1, 'npc'), 
  just = c('right', 'bottom'), name = 'right')
pushViewport(vr)
grid.draw(gt)
upViewport()

dev.off()

##### plot misc GeneFam ortholog map


