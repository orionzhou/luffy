dirI = file.path(DIR_Misc2, "mtgea")
dirO = file.path(DIR_R, "mtgea")
f1 = file.path(dirI, "11_exp.txt")
d1 = read.table(f01, head=TRUE, sep="\t", quote="\"")
f2 = file.path(dirI, "12_cat.txt")
cat1 = read.table(f02, head=TRUE, sep="\t", quote="\"")
cat1[duplicated(cat1$gene),]
d2 = merge(d1, cat1, by.x='gene', by.y='gene', all.x=TRUE)
d3 = d2[, c(1, ncol(d2), c(2:(ncol(d2)-1)))]

e = d3[,-c(1,2)]
filter1 <- function(x, thresh) {sum(x < thresh) == length(x)}
idx = apply(e, 1, filter1, thresh=15)
#idx = apply(e, 1, filter1, thresh=1)

write.table(d3[idx, ], file.path(dirO, "01_filtered_genes.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
d4 = d3[!idx, ]
rownames(d4) = 1:nrow(d4)
e4 = d4[,-c(1,2)]

#select a probe to plot
tissues <- c('leaf', 'petiole', 'stem', 'vegbud', 'flower', 'seed10', 'seed12', 'seed16', 'seed20', 'seed24', 'seed36', 'pod', 'root', 'root0', 'nod4', 'nod10', 'nod14', 'nod28');
tissues.factor <- factor(rep(tissues, each=3), levels=tissues);
plotProbe <- function(exp, tissue, x) { 
  a = data.frame(exp=exp[x,], tissue=tissue);
  p <- ggplot(a, aes(x=tissue, y=exp)) +
    xlab("tissue") + ylab("expression value") +
    geom_boxplot() +
    opts(axis.text.x = theme_text(angle = 45, hjust = 1, size = 6, colour = "grey50"), title=rownames(exp)[x]);
  ggsave(p, filename=paste("probe_", x, ".png", sep=""), path=dirO, width=4, height=3);
}
plotProbe(e4, tissues.factor, 10)

#merge biogical replicates into 1 column (mtGEA data)
colMean <- function(row, num) { apply( matrix(row, num, byrow=FALSE), 2, mean) }
e02 <- t(apply(e01, 1, colMean, num=3))
e02.df <- data.frame(e02)
colnames(e02) = tissues
rownames(e02) = rnames

#for mtGEA fold change analysis
foldchange <- function(exp, group1=c(1:14), group2=c(15:18)) {
  a = max(exp[group2])
  tmp = a >= exp[group1] * 2
  sum(tmp) == length(group1)
}
i10 = apply(e02, 1, foldchange)

ttest_ns <- function(exp, group1, group2) {
  tstat = t.test(exp[group1], exp[group2], alternative="less")
  as.numeric(tstat[1:3]) 
}
t01 <- t(apply(e01, 1, ttest_ns, group1=c(1:42), group2=c(46:54)))
t02 <- t(apply(e01, 1, ttest_ns, group1=c(1:42), group2=c(43:45)))
t03 <- t(apply(e01, 1, ttest_ns, group1=c(1:42), group2=c(46:48)))
t04 <- t(apply(e01, 1, ttest_ns, group1=c(1:42), group2=c(49:51)))
t05 <- t(apply(e01, 1, ttest_ns, group1=c(1:42), group2=c(52:54)))
adjustP <- function(a) {
  b = cbind(a, p.adjust(a[,3], method="bonferroni"), p.adjust(a[,3], method="BH"))
  rownames(b) = rnames
  colnames(b) = c("t", "df", "p_none", "p_bonferroni", "p_fdr")
  b
}
t01 = adjustP(t01)
t02 = adjustP(t02)
t03 = adjustP(t03)
t04 = adjustP(t04)
t05 = adjustP(t05)
t = cbind(t01[,5]<=0.05, t02[,5]<=0.05, t03[,5]<=0.05, t04[,5]<=0.05, t05[,5]<=0.05, i10)
t = cbind(t, (t[,1]|t[,2]|t[,3]|t[,4]|t[,5]) & t[,6])
colnames(t) = c("nod10;nod14;nod28", "nod4", "nod10", "nod14", "nod28", "2-fold", "final")
rownames(t) = rnames
write.table(t+0, file.path(dirI, "21_ttest.txt"), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(rnames[t[,7]], file.path(dirI, "22_nod_up.txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

p <- ggplot(t03) +
  geom_bar(aes(p, count, fill=adjust), position="dodge", stat="identity") +
  opts(axis.text.x=theme_text(size=8, hjust=1, angle=45));
ggsave(p, filename="001.png", path=dirO, width=8, height=6);

#scale data by row for hclust
e5 <- t(scale(t(e4), center=TRUE, scale=TRUE))
d5 <- data.frame(d4[, c(1,2)], e5)

#hclust
e = e5
cordist <- function(x) as.dist(1-cor(t(x), method="pearson"));
e.c.dist <- cordist(t(e));
e.c.hc <- hclust(e.c.dist);
e.c.dendro <- as.dendrogram(e.c.hc);
e.r.dist <- cordist(e);
e.r.hc <- hclust(e.r.dist, method='ward');
e.r.dendro <- as.dendrogram(e.r.hc);

#plot dendrogram and cut trees
hc = e.r.hc
k = 10
png(filename=file.path(dirO, "003.png"), width=800, height=400, units='px');
plot(hc, hang=0.1, cex=0.1, ann=FALSE);
groups <- cutree(hc, k=k);
groups.t <- table(groups)[unique(groups[hc$order])];
m <- c(0, cumsum(groups.t));
for (i in 1:k) {
  rect(m[i] + 0.66, par("usr")[3L], m[i+1]+0.33, mean(rev(hc$height)[(k-1):k]), border=2);
  text(mean(m[i:(i+1)]), mean(c(par("usr")[3L], rev(hc$height)[k])), names(m)[i+1], col='red', cex=1.5);
}
dev.off();

#write dendrogram
dendro = e.r.dendro
tmp = data.frame(treeIdx=1:nrow(e), idx=e.r.hc$order)
sum = data.frame(tmp[order(tmp$idx),], d5[,c(1,2)], subgroup=groups, d5[,-c(1,2)])
sum = sum[, c(2,1,3:ncol(sum))]
write.table(sum, file.path(dirO, "02_summary.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

et = data.frame(as.table(e))
colnames(et) = c("idx", "tissue", "expression")
et2 = merge(et, sum, by.x="idx", by.y="idx")
p = ggplot(et2, aes(treeIdx, tissue)) + geom_tile(aes(fill=expression)) + 
  scale_fill_gradient(low='white', high='red') +
  opts(axis.text.x = theme_text(hjust=1, size=10)) +
  opts(axis.text.y = theme_text(hjust=1, size=10));
ggsave(p, filename = file.path(dirO, "004.png"), width=10, height=5);

f03 <- file.path(DIR_Misc2, 'nbs', '08_anno.txt');
an01 <- read.table(file=f03, sep='\t', quote="\"");
an02 <- unique(as.matrix(an01[,2]));
an02 <- cbind(an02, NA);
for (i in 1:nrow(an02)) {
  probeId = as.character(an02[i, 1]);
  idx = which(a02$id == probeId);
  cat(i, idx, probeId, "\n", sep='\t');
  an02[i, 2] = idx;
}
an02 <- data.frame(an02);
colnames(an02) = c("probeid", "idx");
e05 = e03[an02$idx, ];
rownames(e05) = an02$probeid;

e05.dr <- cordist(e05);
e05.dc <- cordist(t(e05));
e05.fitr <- hclust1(e05.dr);
e05.fitc <- hclust1(e05.dc);
e05.dendro <- as.dendrogram(e05.fitr);
k <- 10;

tree = e05.fitr;
png(filename=file.path(dirO,"008.png"), width=800, height=500, units='px');
plot(tree, hang=0.1, cex=0.7);
cluster <- cutree(tree, k=k);
clustab <- table(cluster)[unique(cluster[tree$order])];
m <- c(0, cumsum(clustab));
for (i in 1:k) {
  rect(m[i] + 0.66, par("usr")[3L], m[i+1]+0.33, mean(rev(tree$height)[(k-1):k]), border=2);
  text(mean(m[i:(i+1)]), mean(c(par("usr")[3L], rev(tree$height)[k])), names(m)[i+1], col='red', cex=1.5);
}
dev.off();

