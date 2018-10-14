require(plyr)
require(ape)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(GenomicRanges)

dirg = file.path(Sys.getenv("genome"), "Zmays_v4")
dirw = file.path(Sys.getenv("misc2"), "grn23")
diro = file.path(Sys.getenv("misc2"), "grn23", "41.qc")

fm = file.path(dirw, "00.0.srr.tsv")
tm = read.table(fm, sep = "\t", header = T, as.is = T)
fg = file.path(dirg, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,16:18)]
gb = group_by(tg, par)
tg2 = summarise(gb, fam = names(sort(table(cat3), decreasing = T))[1])

### hclust
e1 = ti[,-1]
stopifnot(identical(tm$SampleID, colnames(e1)))
colnames(e1) = tm$Tissue

#e5 <- t(scale(t(e4), center=TRUE, scale=TRUE))
#d5 <- data.frame(d4[, c(1,2)], e5)

e = e1
cor_opts = c("spearman", "pearson")
hc_opts = c("ward.D")

cor_opt = "spearman"
hc_opt = "ward.D"
for (cor_opt in cor_opts) {
for (hc_opt in hc_opts) {

plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
#e.r.dist <- cordist(e)
#e.r.hc <- hclust(e.r.dist, method='ward')
#e.r.dendro <- as.dendrogram(e.r.hc)

hc = e.c.hc
fo = sprintf("%s/01.hc.%s.%s.pdf", diro, cor_opt, hc_opt)
pdf(fo, width = 6, height = 9)
#plot(as.dendrogram(e.c.hc, hang = 0.02), cex = 0.6, ann = T, horiz = T)
plot(as.phylo(e.c.hc), cex = 0.8, label.offset = 0.02, no.margin = T)
text(0.001, 65, plot_title, adj = 0)
dev.off()

}
}

### RPM distribution
fr = file.path(dirw, "33.rpm.tsv")
tr = read.table(fr, sep = "\t", header = T, as.is = T)
trl = reshape(tr, direction = 'long', varying = list(2:ncol(tr)), idvar = c("gid"), timevar = "sid", v.names = 'rpm', times = colnames(tr)[2:ncol(tr)])

summary(trl$rpm)
p1 = ggplot(trl) +
  geom_boxplot(aes(x = sid, y = rpm), outlier.shape = NA) + #, draw_quantiles = c(0.25, 0.5, 0.75)) + 
  coord_flip() +
  scale_x_discrete(name = '', breaks = tm$SampleID, labels = tm$Tissue) +
  scale_y_continuous(name = 'RPM', limits = c(0, 30)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))
fp = sprintf("%s/11.rpm.pdf", diro)
ggsave(p1, filename = fp, width = 6, height = 9)

y = apply(tr[,2:ncol(tr)], 2, myfun <- function(x) {
	x1 = sort(x, decreasing = T)
	c('top[01-05]' = sum(x1[1:5]),
		'top[06-10]' = sum(x1[6:10]),
		'top[11-15]' = sum(x1[11:15]),
		'top[16-20]' = sum(x1[16:20]),
		'top[21-30]' = sum(x1[21:30]),
		'top[31-40]' = sum(x1[31:40]),
		'top[41-50]' = sum(x1[41:50])
	)
})
y = cbind(tag = rownames(y), as.data.frame(y))
yl = reshape(y, direction = 'long', varying = list(2:ncol(y)), idvar = c("tag"), timevar = "sid", v.names = 'rpm', times = colnames(y)[2:ncol(y)])

yl$tag = factor(yl$tag, levels = rev(sort(unique(yl$tag))))
p1 = ggplot(yl) +
  geom_bar(aes(x = sid, y = rpm/1000000, fill = tag), position = 'stack', stat = 'identity', width = 0.7) +
  coord_flip() +
  scale_x_discrete(name = '', breaks = tm$SampleID, labels = tm$Tissue) +
  scale_y_continuous(name = 'Proportion Reads', expand = c(0,0), limits = c(0, 1)) +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(1,1,0.1,0.1), "lines")) +
  theme(legend.position = c(0.7, 0.7), legend.direction = "vertical", legend.justification = c(0,0), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))
fp = sprintf("%s/11.rpm.top50.pdf", diro)
ggsave(p1, filename = fp, width = 6, height = 9)

z = ddply(trl, .(sid), myfun <- function(x) {
	x1 = x[order(x[,'rpm'], decreasing = T),][1:10,]
	x2 = cbind.data.frame(x1, crpm = as.numeric(cumsum(x1[,'rpm'])))
	x2
})
z2 = merge(z, tg2, by.x = 'gid', by.y = 'par')
z3 = merge(z2, tm[,c(1,3)], by.x = 'sid', by.y = 'SampleID')
z4 = z3[order(z3$sid, -z3$rpm),]

fo = file.path(diro, '11.rpm.top10.tsv')
write.table(z4, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### RPM correction
yt = ddply(trl, .(sid), summarise, rpm = sum(rpm))
y = apply(tr[,2:ncol(tr)], 2, myfun <- function(x) {
	x1 = sort(x, decreasing = T)
	c(
		'top0' = 0,
		'top5' = sum(x1[1:5]),
		'top10' = sum(x1[1:10]),
		'top20' = sum(x1[1:20]),
		'top30' = sum(x1[1:30]),
		'top50' = sum(x1[1:50]),
		'top100' = sum(x1[1:100])
	)
})
y = cbind.data.frame(opt = rownames(y), y, stringsAsFactors = F)
yl = reshape(y, direction = 'long', varying = list(2:ncol(y)), idvar = c("opt"), timevar = "sid", v.names = 'rpm', times = colnames(y)[2:ncol(y)])


### check housekeeping genes
fh = '/home/springer/zhoux379/data/genome/Zmays_v4/housekeeping/11.tsv'
thk = read.table(fh, sep = "\t", header = T, as.is = T)
hgids = thk$ngid
th = trl[trl$gid %in% hgids,]
th$gid = factor(th$gid, levels = thk$ngid)

opts = c("top0", "top5", "top10", "top20", "top30", "top50", "top100")
opt = opts[7]
z = yl[yl$opt == opt,]
z2 = merge(z, yt, by = 'sid')
z3 = cbind(z2[,c('sid','opt')], sf = z2$rpm.y / (z2$rpm.y - z2$rpm.x))
th2 = merge(th, z3, by = 'sid')
th3 = cbind(th2, rpmc = th2$rpm * th2$sf)

cols = brewer.pal(nrow(thk), 'Set1')
p1 = ggplot(th3) +
  geom_bar(aes(x = sid, y = rpmc, fill = gid), position = 'stack', stat = 'identity', width = 0.7) +
  coord_flip() +
  scale_x_discrete(name = '', breaks = tm$SampleID, labels = tm$Tissue) +
  scale_y_continuous(name = 'RPM', expand = c(0,0)) +
  scale_fill_manual(values = cols, breaks = thk$ngid, labels = thk$name) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,2,0.1,0.1), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0,0), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))
fp = sprintf("%s/13.hk.%s.pdf", diro, opt)
ggsave(p1, filename = fp, width = 7, height = 9)

tb = data.frame()
for (opt in opts) {
	z = yl[yl$opt == opt,]
	z2 = merge(z, yt, by = 'sid')
	z3 = cbind(z2[,c('sid','opt')], sf = z2$rpm.y / (z2$rpm.y - z2$rpm.x))
	th2 = merge(th, z3, by = 'sid')
	th3 = cbind(th2, rpmc = th2$rpm * th2$sf)
	tb = rbind(tb, th3[,c('sid','gid','opt','rpmc')])
}
tb2 = ddply(tb, .(opt, gid), summarize, cv = (sd(rpmc)/mean(rpmc))*100)
tb2$opt = factor(tb2$opt, levels = opts)

cols = brewer.pal(length(opts), 'Paired')
p1 = ggplot(tb2) +
  geom_bar(aes(x = gid, y = cv, fill = opt), position = 'dodge', stat = 'identity', width = 0.8) +
  coord_flip() +
  scale_x_discrete(name = '', breaks = thk$ngid, labels = thk$name) +
  scale_y_continuous(name = 'C.V. of RPM across tissues') +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,2,0.1,0.1), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0,0), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))
fp = sprintf("%s/14.hk.pdf", diro)
ggsave(p1, filename = fp, width = 7, height = 7)


### RPKM distribution
til = reshape(ti, direction = 'long', varying = list(2:ncol(ti)), idvar = c("gid"), timevar = "sid", v.names = 'rpkm', times = colnames(ti)[2:ncol(ti)])

summary(til$rpkm)
p1 = ggplot(til) +
  geom_boxplot(aes(x = sid, y = rpkm), outlier.shape = NA) + #, draw_quantiles = c(0.25, 0.5, 0.75)) + 
  coord_flip() +
  scale_x_discrete(name = '', breaks = tm$SampleID, labels = tm$Tissue) +
  scale_y_continuous(name = 'RPKM', limits = c(0, 30)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))
fp = sprintf("%s/12.rpkm.pdf", diro)
ggsave(p1, filename = fp, width = 6, height = 9)

### Averaging replicates
trl2 = merge(trl, tm[, c("SampleID", "Tissue")], by.x = 'sid', by.y = 'SampleID')
grp = group_by(trl2, Tissue, gid)
trl3 = as.data.frame(summarise(grp, rpm = mean(rpm)))
trw = reshape(trl3, direction = 'wide', timevar = c('Tissue'), idvar = c('gid'))
colnames(trw)[2:ncol(trw)] = gsub("rpm.", "", colnames(trw)[2:ncol(trw)])

fo = file.path(dirw, '35.rpm.tissue.tsv')
write.table(trw, fo, sep = "\t", row.names = F, col.names = T, quote = F)

# RPM stats
ft = file.path(dirw, '35.rpm.tissue.tsv')
tt = read.table(ft, sep = "\t", header = T, as.is = T)
y = apply(tt[,2:ncol(tt)], 2, myfun <- function(x) {quantile(x, seq(0, 1, 0.05)) })
fo = file.path(diro, '15.rpm.tissue.quantile.tsv')
write.table(t(y), fo, sep = "\t", row.names = T, col.names = T, quote = F)

x = apply(tt[,2:ncol(tt)], 2, myfun <- function(x) {sum(x > 0)})
#24,000 - 28,000 

nt = apply(tt[,2:ncol(tt)], 1, myfun <- function(x) {sum(x>0)})
max_expr = apply(tt[,2:ncol(tt)], 1, myfun <- function(x) max(x))
t_nt = table(nt)
tnt = data.frame(ntissue = as.numeric(names(t_nt)), ngene = as.numeric(t_nt))
p = ggplot(tnt) +
  geom_bar(aes(x = ntissue, y = ngene), stat = 'identity', width = 0.6) +
  scale_y_continuous(name = "# genes") +
  scale_x_continuous(name = "# tissues with RPM > 0") + 
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 90, hjust = 0.5))
ggsave(p, filename = file.path(diro, "17.rpm.maf0.pdf"), width = 4, height = 4)

### RPM filtering
sum(nt >= 5)
sum(nt >= 5 & max_expr >= 1)
idxs = which((nt >= 5 & max_expr >= 1) | (nt >= 1 & max_expr >= 3))
length(idxs)
tt2 = tt[idxs,]

fo = file.path(dirw, '36.rpm.filtered.tsv')
write.table(tt2, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### compute RPKM
dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
f_gtb = file.path(dirg, "51.gtb")
f_tbl = file.path(dirg, "51.tbl")
tg = read.table(f_gtb, sep = "\t", header = T, as.is = T)[,1:2]
tt = read.table(f_tbl, sep = "\t", header = F, as.is = T)
colnames(tg) = c("tid", "gid")
colnames(tt) = c("chr", "beg", "end", "srd", "tid", "type", "fam")
tt2 = tt[tt$type %in% c('cds', 'utr5', 'utr3'),]
tg2 = merge(tg, tt2, by = 'tid')

gr = with(tg2, GRanges(seqnames = chr, ranges = IRanges(beg, end = end), gid = gid))
x = unlist(reduce(split(gr, elementMetadata(gr)$gid)))
tr = data.frame(gid = names(x), chr = seqnames(x), beg = start(x), end = end(x), stringsAsFactors = F)
grp = dplyr::group_by(tr, gid)
tr2 = dplyr::summarise(grp, len = sum(end - beg + 1))

fr = file.path(dirw, '36.rpm.filtered.tsv')
to = read.table(fr, sep = "\t", header = T, as.is = T)
to2 = merge(to, tr2, by = 'gid')
stopifnot(nrow(to) == nrow(to2))
for (i in 2:(ncol(to2)-1)) {
	to2[,i] = to2[,i] / (to2$len/1000)
}
to3 = to2[,-ncol(to2)]
dim(to3)
fo = file.path(dirw, '37.rpkm.filtered.tsv')
write.table(to3, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### PCA
fi = file.path(dirw, '36.rpm.filtered.tsv')
fi = file.path(dirw, '37.rpkm.filtered.tsv')
ti = read.table(fi, sep = "\t", header = T, as.is = T)

e1 = ti[,-1]

n_noexp = apply(e1, 1, myfunc <- function(x) sum(x<1))
n_exp = apply(e1, 1, myfunc <- function(x) sum(x>=1))

e = asinh(e1[n_noexp < 10,])
pca <- prcomp(e, center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
tp = cbind.data.frame(SampleID = rownames(x), x[,1:5], stringsAsFactors = F)

p1 = ggplot(tp) +
  geom_point(aes(x = PC1, y = PC2)) +
  geom_text(aes(x = PC1, y = PC2, label = SampleID), nudge_y = 0.01) +
  scale_x_continuous(name = 'PC1') +
  scale_y_continuous(name = 'PC2') +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fp = sprintf("%s/01.sample.pca.pdf", diro)
ggsave(p1, filename = fp, width = 7, height = 7)

### LDA
require(MASS)

fm = file.path(dirw, "00.0.srr.tsv")
tm = read.table(fm, sep = "\t", header = T, as.is = T)
fr = file.path(dirw, "33.rpm.tsv")
tr = read.table(fr, sep = "\t", header = T, as.is = T)
tissue_map = tm$Tissue; names(tissue_map) = tm$SampleID

e1 = tr[,-1]
n_noexp = apply(e1, 1, myfunc <- function(x) sum(x<1))
idxs = which(n_noexp < 10)
length(idxs)

e = asinh(e1[n_noexp < 10,])
pca <- prcomp(e, center = F, scale. = F)

prop.pca = pca$sdev^2/sum(pca$sdev^2)
prop.pca[1:5]

pca.x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
tp = cbind.data.frame(tissue = tissue_map[rownames(x)], 
	pca.x = x[,1], pca.y = x[,2], stringsAsFactors = F)

p2 <- ggplot(tp) + geom_point(aes(pca.x, pca.y, colour = tissue), size = 2) +
  labs(x = paste("PC1 (", prop.pca[1]*100, ")", sep=""),
       y = paste("PC2 (", prop.pca[2]*100, ")", sep=""))
fp = sprintf("%s/01.sample.pca.pdf", diro)
ggsave(p2, filename = fp, width = 10, height = 7)


tl = asinh(t(tr[idxs,-1]))
colnames(tl) = tr$gid[idxs]
tl = data.frame(SampleID = rownames(tl), tl, stringsAsFactors = F)
tl2 = merge(tm[,c('SampleID', 'Tissue')], tl, by = 'SampleID')
tl2[1:5,1:5]
tl3 = tl2[,-1]


r <- lda(formula = Tissue ~ ., data = tl3)#, prior = c(1,1,1)/3)
r$svd^2/sum(r$svd^2)

r$prior
r$counts
#r$means
#r$scaling
prop.lda = r$svd^2/sum(r$svd^2)
prop.lda

plda <- predict(object = r, newdata = tl3)
dataset = data.frame(tissue = tl3$Tissue, lda = plda$x)
                     
p1 <- ggplot(dataset) + 
	geom_point(aes(lda.LD1, lda.LD2, colour = tissue), size = 2) + 
  	labs(x = paste("LD1 (", prop.lda[1]*100, ")", sep=""),
       y = paste("LD2 (", prop.lda[2]*100, ")", sep=""))
fp = sprintf("%s/01.sample.lda.pdf", diro)
ggsave(p1, filename = fp, width = 10, height = 7)
