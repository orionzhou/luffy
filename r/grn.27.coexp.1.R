source("br.fun.R")

#dirw = file.path(Sys.getenv("misc2"), "grn23", "47.coexp.test")
dirw = file.path(Sys.getenv("misc2"), "briggs", "47.coexp.test")

#fi = file.path(dirw, "../37.rpkm.filtered.tsv")
fi = file.path(dirw, "../36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)
ti = spread(ti[ti$Genotype == "B73", c(1,2,5)], Tissue, fpkm)

expr = t(as.matrix(ti[,-1]))
colnames(expr) = ti[,1]
gids = ti[,1]
ng = length(gids)


opt = 3
## camoco
if (opt == 1) {

pcc.matrix = cor(asinh(expr), method = 'pearson')
pcc = pcc.matrix[lower.tri(pcc.matrix)]

pcc[pcc == 1] = 0.999999
pcc[pcc == -1] = -0.999999
#pcc2 = log((1+pcc) / (1-pcc)) / 2
pcc2 = atanh(pcc)
coexv = (pcc2 - mean(pcc2)) / sd(pcc2)

coexm <- matrix(rep(1, ng*ng), nrow=ng)
coexm[lower.tri(coexm)] = coexv
coexm = t(coexm)
coexm[lower.tri(coexm)] = coexv

fp = sprintf("%s/01.edgeweight/C.rda", dirw)
save(coexv, coexm, file = fp)

fp = sprintf("%s/01.edgeweight/C.pdf", dirw)
pdf(fp, width = 8, height = 8)
hist(coexv, 100)
dev.off()

}

## clusterOne
opt = 3
if (opt == 2) {

pcc.matrix = cor(asinh(expr), method = 'pearson')
pcc = pcc.matrix[lower.tri(pcc.matrix)]

coex <- matrix(rep(-1, ng*ng), nrow=ng)
coex[lower.tri(coex)] = pcc
coex = t(coex)
coex[lower.tri(coex)] = pcc

dummy <- matrix(rep(NA, ng*ng), nrow=ng)
ii = which(lower.tri(dummy))
colidx = as.integer((ii - 1) / ng) + 1
rowidx = ii - (colidx - 1) * ng

rankmatrix = t(apply(-coex, 1, rank))
coexv = sapply(1:length(pcc), myfunc <- function(i) 
	sqrt(rankmatrix[rowidx[i], colidx[i]] * rankmatrix[colidx[i], rowidx[i]])
)

coexm <- matrix(rep(1, ng*ng), nrow=ng)
coexm[lower.tri(coexm)] = coexv
coexm = t(coexm)
coexm[lower.tri(coexm)] = coexv

fp = sprintf("%s/01.edgeweight/R.rda", dirw)
save(coexv, coexm, file = fp)

fp = sprintf("%s/01.edgeweight/R.pdf", dirw)
pdf(fp, width = 8, height = 8)
hist(coexv, 100)
dev.off()

}

## WGCNA
opt = 3
if (opt == 30) {

powers = c(c(1:10), seq(from = 12, to=20, by=2))
powers = 1:20
sft = pickSoftThreshold(expr, powerVector = powers, verbose = 5)
fo = sprintf("%s/01.edgeweight/01.power.pdf", dirw)
pdf(fo, width = 9, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

}

if (opt == 3) {

enableWGCNAThreads()
allowWGCNAThreads()

softPower = 8#10
adjacency = adjacency(expr, power = softPower, type = "signed", corFnc = "cor")
dim(adjacency)

TOM = TOMsimilarity(adjacency, TOMType = "signed")
coexm = TOM
coexv = coexm[lower.tri(coexm)]

fp = sprintf("%s/01.edgeweight/W.rda", dirw)
save(coexv, coexm, file = fp)

fp = sprintf("%s/01.edgeweight/W.pdf", dirw)
pdf(fp, width = 8, height = 8)
hist(coexv, 100)
dev.off()

}