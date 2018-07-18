source("br.fun.R")
require(Hmisc)
require(tidyr)
require(RColorBrewer)

dirw = "/home/springer/zhoux379/data/misc1/li2013"
chrs = sprintf("%d", 1:10)

### DE of trans-eQTL hotspots
ff = '/home/springer/zhoux379/data/misc2/briggs/43.deseq/01.de.tsv'
tf = read.table(ff, header = T, sep = "\t", as.is = T)
tf = tf[tf$comp == 'B73 vs Mo17', -2]

tk = tf
tk$is.de[tk$is.de=='up'] = 1
tk$is.de[tk$is.de=='down'] = -1
tk$is.de = as.integer(tk$is.de)
tk = spread(tk, tissue, is.de)
tk[is.na(tk)] = 0


fh = file.path(dirw, '49.hs.module.tsv')
th = read.table(fh, header = T, sep = "\t", as.is = T, quote = '')
th$opt = sapply(strsplit(th$mid, split = "[.]"), '[', 1)

th2 = merge(th[,1:2], tk, by='gid')
th2 = within(th2, {gid = sprintf("%s.%s", mid, gid)})

th3 = data.frame()
for (hs in sort(unique(th2$mid))) {
	ths = th2[th2$mid == hs,]
	hcl = hclust(dist(ths[,-c(1:2)]), method = "ward.D")
	th3 = rbind(th3, ths[hcl$order,])
}

th3 = cbind(y = 1:nrow(th3), th3)

ty = ddply(th3, .(mid), summarise, ymin=min(y), ymax=max(y), pos=(min(y)+max(y))/2)

tp = gather(th3[,-c(2:3)], tissue, is.de, -y)
tp$is.de = as.character(tp$is.de)
tp$tissue = factor(tp$tissue, levels = tissues)

p1 = ggplot(tp) +
  geom_tile(aes(x = tissue, y = y, fill = is.de)) + 
  geom_segment(data = ty, aes(x=0.3,xend=0.5,y=ymin,yend=ymin)) +
  geom_segment(data = ty, aes(x=0.3,xend=0.5,y=ymax,yend=ymax)) +
  scale_x_discrete(name = '', expand=c(0,0)) +
  scale_y_continuous(name = 'Genes', expand= c(0,0), breaks=ty$pos, labels=ty$mid) +
  scale_fill_manual(values = c("white", "#E41A1C", "#4DAF4A"), labels = c("Not DE", "B < M", "B > M")) + 
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill='grey', size=1)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 70, hjust = 1)) +
  theme(axis.text.y = element_text(size=8))
fo = sprintf("%s/31.de.pdf", dirw)
ggsave(p1, filename = fo, width = 5, height = 11)

## module eigen-gene and eigen-log2(B/H)
dirx = '/home/springer/zhoux379/data/misc2/briggs'
ff = file.path(dirx, '36.long.filtered.tsv')
tf = read.table(ff, header = T, sep = "\t", as.is = T)
tf = tf[tf$Genotype != 'Mo17xB73',]
tf$Genotype = factor(tf$Genotype, levels = c("B73","Mo17", "B73xMo17"))

fh = file.path(dirw, '49.hs.module.tsv')
th = read.table(fh, header = T, sep = "\t", as.is = T, quote = '')
th = th[th$gid %in% unique(tf$gid),]
th$opt = sapply(strsplit(th$mid, split = "[.]"), '[', 1)

### look at ME
tme = data.frame()
for (gt in gts) {
	tiw = spread(tf[tf$Genotype == gt, c('gid','Tissue','fpkm')], Tissue, fpkm)
	tw = merge(th[,1:2], tiw, by = 'gid')
	cols = tw$mid

	datExpr = t(as.matrix(tw[,-c(1:2)]))
	colnames(datExpr) = sprintf("%s.%s", tw$mid, tw$gid)
	me = moduleEigengenes(datExpr, cols)$eigengenes
	tme1 = data.frame(t(me))
	colnames(tme1) = colnames(tiw)[-1]
	tme1 = cbind(gt = gt, mid = colnames(me), tme1)
	tme = rbind(tme, tme1)
}
tme$mid = substr(tme$mid, 3, nchar(tme$mid))
opts = sapply(strsplit(tme$mid, split = "[.]"), '[', 1)

mids_order = c()
for (opt in sort(unique(opts))) {
    tme1 = tme[opts == opt & tme$gt == 'B73',]
	hcl = hclust(dist(tme1[,-c(1:3)]), method = "ward.D")
	mids_order = c(mids_order, tme1$mid[hcl$order])
}
tmid = data.frame(mid = mids_order, y = 1:length(mids_order), stringsAsFactors = F)
tmid2 = cbind(tmid, opt = sapply(strsplit(tmid$mid, split = "[.]"), '[', 1))
ty = ddply(tmid2, .(opt), summarise, ymin = min(y), ymax = max(y), ypos = (min(y)+max(y))/2)
ty = cbind(ty, lab = paste(substr(ty$opt, nchar(ty$opt), nchar(ty$opt)), " %up% ."))

tme = merge(tme, tmid, by = 'mid')

tml = gather(tme, tissue, me, -gt, -mid, -y)
tml$tissue = factor(tml$tissue, levels = tissues)
tml$gt = factor(tml$gt, levels = gts)

p1 = ggplot(tml) +
  geom_tile(aes(x = tissue, y = y, fill = me)) + 
  #geom_segment(data = ty, aes(x=0.0,xend=0.5,y=ymin,yend=ymin)) +
  #geom_segment(data = ty, aes(x=0.0,xend=0.5,y=ymax,yend=ymax)) +
  geom_segment(data = ty, aes(x=0.5,xend=0.5,y=ymin,yend=ymax)) +
  scale_x_discrete(name = '', expand=c(0,0)) +
  scale_y_continuous(name = 'trans-eQTL hotspot', expand=c(0,0), breaks=ty$ypos, labels=parse(text = ty$lab)) +
  #scale_fill_manual(values = ) + 
  scale_fill_gradient2() +
  facet_wrap(~gt) + 
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill='grey', size=1)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 70, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8))
fo = sprintf("%s/32.module.eigen.pdf", dirw)
ggsave(p1, filename = fo, width = 10, height = 6)

### look at log2(B/M)
th2 = merge(th[,1:2], tf, by = 'gid')
th3 = spread(th2[,-6], Genotype, fpm)
th3 = within(th3, {bmr = log(B73/Mo17)})
describe(th3$bmr)
th3$bmr[th3$bmr == Inf] = 6
th3$bmr[th3$bmr == -Inf] = -6
th3$bmr[is.nan(th3$bmr)] = 0

th4 = spread(th3[,-c(4:6)], Tissue, bmr)

me = moduleEigengenes(t(th4[,-c(1:2)]), th4$mid)$eigengenes
tme = data.frame(t(me))
colnames(tme) = colnames(th4)[-c(1:2)]
rownames(tme) = substr(rownames(tme), 3, nchar(rownames(tme)))

mids = sort(unique(th4$mid))
to = data.frame(matrix(NA, nrow=length(mids), ncol=length(tissues)))
rownames(to) = mids
colnames(to) = colnames(th4)[-c(1:2)]
for (mid in mids) {
    tw1 = th4[th4$mid == mid,]
    e = tw1[,-c(1:2)]
    e = t(scale(t(e)))
    svd1 = svd(e, nu = 1, nv = 1)
    pc1 = svd1$v[,1]

    scaledExpr = scale(t(e))
    averExpr = rowMeans(scaledExpr, na.rm = TRUE)
    corAve = cor(averExpr, pc1, use = "p")
    if (!is.finite(corAve)) corAve = 0;
    if (corAve<0) pc1 = -pc1
    stopifnot(sum(pc1-tme[mid,]) == 0)
    to[mid,] = pc1
}
to = cbind(mid = rownames(to), to)
opts = sapply(strsplit(to$mid, split = "[.]"), '[', 1)

mids_order = c()
for (opt in sort(unique(opts))) {
    to1 = to[opts == opt,]
	hcl = hclust(dist(to1[,-c(1:3)]), method = "ward.D")
	mids_order = c(mids_order, to1$mid[hcl$order])
}
tmid = data.frame(mid = mids_order, y = 1:length(mids_order), stringsAsFactors = F)
tmid2 = cbind(tmid, opt = sapply(strsplit(tmid$mid, split = "[.]"), '[', 1))
ty = ddply(tmid2, .(opt), summarise, ymin = min(y), ymax = max(y), ypos = (min(y)+max(y))/2)
ty = cbind(ty, lab = paste(substr(ty$opt, nchar(ty$opt), nchar(ty$opt)), " %up% ."))

to = merge(to, tmid, by = 'mid')
tol = gather(to, Tissue, bpr, -mid, -y)
tol$Tissue = factor(tol$Tissue, levels = tissues)

p1 = ggplot(tol) +
  geom_tile(aes(x = Tissue, y = y, fill = bpr)) + 
  geom_segment(data = ty, aes(x=0.5,xend=0.5,y=ymin,yend=ymax)) +
  scale_x_discrete(name = '', expand=c(0,0)) +
  scale_y_continuous(name = 'trans-eQTL hotspot', expand=c(0,0), breaks=ty$ypos, labels=parse(text = ty$lab)) +
  #scale_fill_manual(values = ) + 
  scale_fill_gradient2(name='log(B/M)', breaks=c(-0.5,0.5), labels=c("B < M", "B > M")) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_text(size=9), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size=8), legend.background = element_rect(fill='grey', size=1)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 70, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8))
fo = sprintf("%s/33.bpr.eigen.pdf", dirw)
ggsave(p1, filename = fo, width = 4, height = 6)

# look at DOA
th2 = merge(th[,1:2], tf, by = 'gid')
th3 = spread(th2[,-6], Genotype, fpm)
th3 = within(th3, {doa = log((B73xMo17 - (B73+Mo17)/2) / (pmax(B73,Mo17)-(B73+Mo17)/2))})
describe(th3$doa)
th3$doa[th3$doa == Inf] = 7
th3$doa[th3$doa == -Inf] = -7
th3$doa[is.nan(th3$doa)] = 0

th4 = spread(th3[,-c(4:6)], Tissue, doa)

me = moduleEigengenes(t(th4[,-c(1:2)]), th4$mid)$eigengenes
tme = data.frame(t(me))
colnames(tme) = colnames(th4)[-c(1:2)]
rownames(tme) = substr(rownames(tme), 3, nchar(rownames(tme)))

to = cbind(mid = rownames(tme), tme)
opts = sapply(strsplit(to$mid, split = "[.]"), '[', 1)

mids_order = c()
for (opt in sort(unique(opts))) {
    to1 = to[opts == opt,]
	hcl = hclust(dist(to1[,-c(1:3)]), method = "ward.D")
	mids_order = c(mids_order, to1$mid[hcl$order])
}
tmid = data.frame(mid = mids_order, y = 1:length(mids_order), stringsAsFactors = F)
tmid2 = cbind(tmid, opt = sapply(strsplit(tmid$mid, split = "[.]"), '[', 1))
ty = ddply(tmid2, .(opt), summarise, ymin = min(y), ymax = max(y), ypos = (min(y)+max(y))/2)
ty = cbind(ty, lab = paste(substr(ty$opt, nchar(ty$opt), nchar(ty$opt)), " %up% ."))

to = merge(to, tmid, by = 'mid')
tol = gather(to, Tissue, bpr, -mid, -y)
tol$Tissue = factor(tol$Tissue, levels = tissues)

p1 = ggplot(tol) +
  geom_tile(aes(x = Tissue, y = y, fill = bpr)) + 
  geom_segment(data = ty, aes(x=0.5,xend=0.5,y=ymin,yend=ymax)) +
  scale_x_discrete(name = '', expand=c(0,0)) +
  scale_y_continuous(name = 'trans-eQTL hotspot', expand=c(0,0), breaks=ty$ypos, labels=parse(text = ty$lab)) +
  #scale_fill_manual(values = ) + 
  scale_fill_gradient2(name='log(DOA)', breaks=c(-0.5,0.5), labels=c("F1 < MP", "F1 > MP")) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_text(size=9), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size=8), legend.background = element_rect(fill='grey', size=1)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 70, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8))
fo = sprintf("%s/34.doa.eigen.pdf", dirw)
ggsave(p1, filename = fo, width = 4, height = 6)
