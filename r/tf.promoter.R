require(dplyr)
require(GenomicRanges)
require(ggplot2)
require(tidyr)

diri = "/home/springer/zhoux379/data/misc2/briggs"
dirw = "/home/springer/zhoux379/data/misc1/tf.promoter"

### TF expression
fx = file.path(diri, '35.long.tsv')
tx = read.table(fx, header = T, sep = "\t", as.is = T)

fi = file.path(dirw, 'tf_v4.tsv')
ti = read.table(fi, header = F, sep = "\t", as.is = T)

gids = ti$V2
gids = unique(gids[gids != ''])

tx2 = tx[tx$gid %in% gids,]

tx3 = cbind(tx2, treat = sprintf("%s|%s", tx2$Genotype, tx2$Tissue))
to1 = spread(tx3[,c('gid','treat','fpm')], treat, fpm)
to2 = spread(tx3[,c('gid','treat','fpkm')], treat, fpkm)
fo = file.path(dirw, "tf_v4_fpm.tsv")
write.table(to1, fo, sep = "\t", row.names = F, col.names = T, quote = F)

fo = file.path(dirw, "tf_v4_fpkm.tsv")
write.table(to2, fo, sep = "\t", row.names = F, col.names = T, quote = F)


### promoter expression
fx = file.path(diri, '35.long.tsv')
tx = read.table(fx, header = T, sep = "\t", as.is = T)

fi = file.path(dirw, '11.promoter.tsv')
ti = read.table(fi, header = F, sep = "\t", as.is = T)

gids = ti$V2
gids = unique(gids[gids != ''])
nrow(ti)
length(gids)

tx2 = tx[tx$gid %in% gids,]

tx3 = cbind(tx2, treat = sprintf("%s|%s", tx2$Genotype, tx2$Tissue))
to1 = spread(tx3[,c('gid','treat','fpm')], treat, fpm)
to2 = spread(tx3[,c('gid','treat','fpkm')], treat, fpkm)
fo = file.path(dirw, "12.promoter.fpm.tsv")
write.table(to1, fo, sep = "\t", row.names = F, col.names = T, quote = F)

fo = file.path(dirw, "12.promoter.fpkm.tsv")
write.table(to2, fo, sep = "\t", row.names = F, col.names = T, quote = F)


## get gene IDs for TFs and promoters
fi = file.path(dirw, "tf_v4_fpkm.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
gids1 = ti$gid

fi = file.path(dirw, "12.promoter.fpkm.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
gids2 = ti$gid

## read in FPKM values for all filtered genes - note that genes expressed in any of the tissues in B73 or Mo17 are filtered and not in this list
fi = file.path(diri, "36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

# check if any TFs or promoters were filtered due to non-expression
length(gids1)
gids1 = gids1[gids1 %in% ti$gid]
length(gids1)

length(gids2)
gids2 = gids2[gids2 %in% ti$gid]
length(gids2)

# compute and normalize co-expression network for all genes
	tiw = spread(ti[ti$Genotype == 'B73', -c(3,4)], Tissue, fpkm)
	expr = t(as.matrix(tiw[,-1]))
	colnames(expr) = tiw[,1]
	gids = tiw[,1]
	ng = length(gids)

	pcc.matrix = cor(asinh(expr), method = 'pearson')
	pcc = pcc.matrix[lower.tri(pcc.matrix)]

	ii = which(lower.tri(pcc.matrix))
	colidx = as.integer((ii - 1) / nrow(pcc.matrix)) + 1
	rowidx = ii - (colidx - 1) * nrow(pcc.matrix)

	pcc[pcc == 1] = 0.999999
	pcc[pcc == -1] = -0.999999
	#pcc2 = log((1+pcc) / (1-pcc)) / 2
	pcc2 = atanh(pcc)
	coexv = (pcc2 - mean(pcc2)) / sd(pcc2)

# extract co-expression network for selected genes (TFs and promoters)
idxs = which(gids %in% c(gids1, gids2))
idxsv = which(colidx %in% idxs & rowidx %in% idxs)

to = data.frame(gid1 = gids[colidx[idxsv]], gid2 = gids[rowidx[idxsv]], z.pcc = coexv[idxsv], stringsAsFactors = F)
fo = file.path(dirw, "31.pcc.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)


### visualize co-expression TF&promoter co-expression network
# for more check: http://kateto.net/network-visualization
require(igraph)
require(network)
#require(sna)

fc = file.path(dirw, "31.pcc.tsv")
tc = read.table(fc, header = T, sep = "\t", as.is = T)

net = graph_from_data_frame(
	d = tc[tc$z.pcc >= 3,], 
	vertices = unique(c(tc$gid1, tc$gid2)),
	directed = F
)
V(net)$color = rep('tomato', length(V(net)))
V(net)$color[which(V(net)$name %in% gids2)] = rep('royalblue', length(gids2))

coords = layout_in_circle(net, order = c(gids1, gids2))

fp = file.path(dirw, "35.network3.pdf")
pdf(fp, width = 8, height = 8)
par(mar=c(8,6,6,6))
plot(net, vertex.size = 3, vertex.label = NA, layout = coords)

x = coords[,1]*1.21
y = coords[,2]*1.21

angle = ifelse(atan(-(x/y))*(180/pi) < 0,  90 + atan(-(x/y))*(180/pi), 270 + atan(-x/y)*(180/pi))

for (i in 1:length(x)) {
	text(x=x[i], y=y[i], labels=V(net)$name[i], adj=NULL, pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
}
dev.off()
