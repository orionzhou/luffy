require(plyr)
require(ape)
require(ggplot2)

options(stringsAsFactors = FALSE)

dirw = file.path(Sys.getenv("misc2"), "grn23", "51.camoco")

dirg = file.path(Sys.getenv("genome"), "Zmays_v4")
fg = file.path(dirg, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,16:18)]
gb = group_by(tg, par)
tg2 = summarise(gb, fam = names(sort(table(cat3), decreasing = T))[1])

### run camoco on a lab workstation
source activate camoco
camoco --help
camoco build-refgen $genome/Zmays_v4/51.gff maize v34 v34 maize
camoco build-cob --rawtype RNASEQ --index-col 1 --min-single-sample-expr 0 $misc2/grn23/50.camoco.tsv grn23 1.0 maize
camoco health --out $misc2/grn23/51.camoco.health/grn23 grn23

source activate camoco
ipython
import camoco as co
x = co.COB("grn23")

x.clusters.to_csv("/home/springer/zhoux379/data/misc2/grn23/51.camoco/clusters.csv")
x.coex.score.to_csv("/home/springer/zhoux379/data/misc2/grn23/51.camoco/coex.csv", index = False, float_format='%g')

### run grn.4.hpc.R

#####
fx = file.path(dirw, "camoco.rda")
x = load(fx)
x



### obtain camoco PCC from scratch
fi = file.path(dirw, "../37.rpkm.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

expr = t(as.matrix(ti[,-1]))
colnames(expr) = ti[,1]

pcc.matrix = cor(asinh(expr), method = 'pearson')
pcc = pcc.matrix[lower.tri(pcc.matrix)]
pcc[pcc == 1] = 0.9999
pcc[pcc == -1] = -0.9999
#pcc2 = log((1+pcc) / (1-pcc)) / 2
pcc2 = atanh(pcc)
pcc3 = (pcc2 - mean(pcc2)) / sd(pcc2)
head(pcc3)

ng = nrow(ti)
coex <- matrix(rep(1, ng*ng), nrow=ng)
coex[lower.tri(coex)] = pcc3
coex = t(coex)
coex[lower.tri(coex)] = pcc3
save(coex, file = file.path(dirw, "11.camoco.dudu.rda"))


### compare camoco and wgcna
dirw = file.path(Sys.getenv("misc2"), "grn23")
fi = file.path(dirw, "37.rpkm.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

fc1 = file.path(dirw, "51.camoco/11.camoco.dudu.rda")
x = load(fc1)
coex1 = coex

fc2 = file.path(dirw, "52.wgcna/x.RData")
x = load(fc2)
coex2 = TOM
tom = coex2[low.tri(coex2)]
pdf(file.path(dirw, "52.wgcna/coex.pdf"), width = 6, height = 6)
hist(tom, breaks = seq(0, 1, by = 0.025))
dev.off()

ec = sapply(1:nrow(coex1), myfunc <- function(i) cor(coex1[i, -i], coex2[i, -i]))
idxs = c(1:length(ec))[order(ec, decreasing = T)][1:100]
novlps = sapply(idxs, myfunc <- function(i) {
	idxs1 = which(coex1[i,-i] >= sort(coex1[i,-i], decreasing = T)[1000])
	idxs2 = which(coex2[i,-i] >= sort(coex2[i,-i], decreasing = T)[1000])
	sum(idxs1 %in% idxs2)
})
novlps = sapply(1:length(ec), myfunc <- function(i) {
	idxs1 = which(coex1[i,-i] >= sort(coex1[i,-i], decreasing = T)[1000])
	idxs2 = which(coex2[i,-i] >= sort(coex2[i,-i], decreasing = T)[1000])
	sum(idxs1 %in% idxs2)
})
cor.test(ec, novlps)

tp = data.frame(gid = ti$gid, ec = ec, novlp = novlps, stringsAsFactors = F)
tp = merge(tp, tg2, by.x = 'gid', by.y = 'par')
tp = tp[order(tp$ec),]
save(tp, file = file.path(dirw, "22.rda"))

p1 = ggplot(tp) +
  geom_point(aes(x = ec, y = novlp), size = 0.5) +
  scale_x_continuous(name = 'Expression Conservation (EC) Score') +
  scale_y_continuous(name = '# Common Neighbors btw. Networks', limits = c(0, 1000)) + 
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, color = "blue")) +
  theme(axis.text.y = element_text(size = 8, color = "blue", angle = 90, hjust  = 0.5))
fp = file.path(dirw, "51.camoco/21.ec.pdf")
ggsave(p1, filename = fp, width = 6, height = 6)

x = load(file.path(dirw, "22.rda"))
require(MASS)
p2 = ggplot(tp, aes(x = ec, y = novlp)) +
  stat_density2d(aes(fill = ..density..), contour = F, geom="tile") +
  scale_x_continuous(name = 'Expression Conservation (EC) Score', limits = c(-1, 1), expand = c(0, 0)) +
  scale_y_continuous(name = '# Common Neighbors btw. Networks', limits = c(0, 1000), expand = c(0, 0)) + 
  scale_fill_gradient(low = "white", high = "royalblue") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.85), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA), legend.key.size = unit(1, 'lines'), legend.title = element_text(size = 8, angle = 0), legend.text = element_text(size = 8, angle = 0)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, color = "blue")) +
  theme(axis.text.y = element_text(size = 8, color = "blue", angle = 90, hjust  = 0.5))
fp = file.path(dirw, "51.camoco/22.ec.density.pdf")
ggsave(p2, filename = fp, width = 6, height = 6)
