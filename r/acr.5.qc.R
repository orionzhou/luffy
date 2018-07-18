require(tidyverse)
require(Hmisc)
require(ape)
require(gtable)
require(grid)
require(RColorBrewer)

dirw = '/home/springer/zhoux379/data/misc1/maize.acr'

### process expression matrix
fi = '/home/springer/nosha003/wgbs_schmitz/ACR/rpm_matrix.txt'
ti = read.table(fi, sep = "\t", header = T, stringsAsFactors = F)
colnames(ti)[1] = 'gid'
ti$gid = sapply(strsplit(ti$gid, split = "[:]"), "[", 2)

apply(ti[,-1], 2, sum)
#for (i in 2:ncol(ti)) {
#	ti[,i] = ti[,i] * 1000000 / sum(ti[,i])
#}
#apply(ti[,-1], 2, sum)

tl = gather(ti, cond, rpm, -gid)
tmp = strsplit(tl$cond, split = "_")
tl = data.frame(tl[,-2], genotype = sapply(tmp, "[", 1),
	tissue = sapply(tmp, "[", 2),
	rep = sapply(tmp, "[", 3), stringsAsFactors = F)
tl$genotype[tl$genotype=='X207'] = 'PH207'

fo = file.path(dirw, "10.rpm.tsv")
write.table(tl, fo, sep = "\t", row.names = F, col.names = T, quote = F)


### read in RPM matrix
fi = file.path(dirw, '10.rpm.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

ti = cbind(ti, cond = sprintf("%s_%s_%s", ti$genotype, ti$tissue, ti$rep))
tw = spread(ti[,-c(3:5)], cond, rpm)

cnames = colnames(tw)[-1]
tmp = strsplit(cnames, split = "_")
tm = data.frame(cond = cnames, genotype = sapply(tmp, "[", 1), tissue = sapply(tmp, "[", 2), rep = sapply(tmp, "[", 3), stringsAsFactors = F)

e1 = tw[,-1]
stopifnot(identical(tm$cond, colnames(e1)))
dim(e1)

n_exp = apply(e1, 1, myfunc <- function(x) sum(x>=1))
e = e1[n_exp >= ncol(e1) * 0.8,]
dim(e)


# hc 
cor_opt = "spearman"
hc_opt = "ward.D"

plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
#e.r.dist <- cordist(e)
#e.r.hc <- hclust(e.r.dist, method='ward')
#e.r.dendro <- as.dendrogram(e.r.hc)

hc = e.c.hc
fo = sprintf("%s/12.hc.%s.%s.pdf", dirw, cor_opt, hc_opt)
pdf(fo, width = 4, height = 4)
#plot(as.dendrogram(e.c.hc, hang = 0.02), cex = 0.6, ann = T, horiz = T)
plot(as.phylo(e.c.hc), cex = 0.8, label.offset = 0.02, no.margin = T)
#text(0.001, 125, plot_title, adj = 0)
dev.off()


### PCA
#e = asinh(e)
pca <- prcomp(e, center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
tp = cbind.data.frame(cond = rownames(x), x[,1:5], stringsAsFactors = F)
tp2 = merge(tp, tm, by = 'cond')
#tp2$tissue = factor(tp2$tissue, levels = unique(tp2$tissue))
#tp2$genotype = factor(tp2$genotype, levels = unique(tp2$genotype))

cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

p1 = ggplot(tp2) +
  geom_point(aes(x = PC1, y = PC2, shape = genotype, color = tissue), size = 3) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(name = ylab) +
  scale_color_manual(name = "", values = cols) +
  scale_shape_manual(name = "", values = c(1,0,2,6)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fp = sprintf("%s/11.sample.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 5, height = 4)



### FPM distribution
fi = file.path(dirw, '10.rpm.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
gts = c("B73", "PH207")
ti = ti[ti$genotype %in% gts,]

# density distribution
to1 = ti
to1 = cbind(to1, asinh_rpm = asinh(to1$rpm))
#describe(to1$asinh_fpm)
to1 = to1[to1$rpm > 0, ]
to1$asinh_rpm[to1$asinh_rpm > 10] = 10

to2 = within(to1, {
    rpm_bin = cut(asinh_rpm, breaks=seq(0,10,0.1), include.lowest = T, labels = seq(0.05, 9.95, by = 0.1))
})
sum(is.na(to2$rpm_bin))
to2$rpm_bin = as.numeric(as.character(to2$rpm_bin))
grp = group_by(to2, tissue, genotype, rpm_bin)
to3 = summarise(grp, num_genes = n(), total_rpm = sum(rpm))
to4 = gather(to3, type, value, -tissue, -genotype, -rpm_bin)

labs = c(0, 1, 10, 100, 1000, 10000)
brks = asinh(labs)

p1 = ggplot(to3) +
  geom_bar(aes(x = rpm_bin, y = num_genes, fill = genotype), stat = 'identity', position = 'dodge') +
  #geom_line(aes(color = genotype, linetype = genotype)) +
  #geom_point(aes(color = genotype, shape = genotype), size = 0.5) +
  scale_x_continuous(name = 'RPM', limits = c(0, 10), breaks = brks, labels = labs) +
  #scale_y_continuous(name = '#') +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap( ~tissue, ncol = 1, scale = 'free') +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.5,0), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, color = "black", angle = 0, hjust = 1))

p2 = ggplot(to3) +
  geom_bar(aes(x = rpm_bin, y = total_rpm, fill = genotype), stat = 'identity', position = 'dodge') +
  #geom_line(aes(color = genotype, linetype = genotype)) +
  #geom_point(aes(color = genotype, shape = genotype), size = 0.5) +
  scale_x_continuous(name = 'RPM', limits = c(0, 10), breaks = brks, labels = labs) +
  #scale_y_continuous(name = '#') +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap( ~tissue, ncol = 1, scale = 'free') +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.5,0), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, color = "black", angle = 0, hjust = 1))

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
gs = list(g1, g2)
wids = c(4,4)
g <- gtable_matrix(name = 'demo', grobs = matrix(gs, ncol = length(gs)), widths = wids, heights = 1)

fp = sprintf("%s/14.rpm.pdf", dirw)
pdf(file = fp, width = 8, height = 6, bg = 'transparent')
grid.newpage()
grid.draw(g)
dev.off()


