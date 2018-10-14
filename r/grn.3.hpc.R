require(plyr)
require(ape)
require(ggplot2)
require(WGCNA)

options(stringsAsFactors = FALSE)

dirw = file.path(Sys.getenv("misc2"), "grn23")
diro = file.path(Sys.getenv("misc2"), "grn23", "52.wgcna")

fi = file.path(dirw, "37.rpkm.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

enableWGCNAThreads()
allowWGCNAThreads()

datExpr = t(as.matrix(ti[,-1]))
colnames(datExpr) = ti[,1]

### hclust and cut tree into modules
softPower = 10
adjacency = adjacency(datExpr, power = softPower)
dim(adjacency)

TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")

#as.dist(1-cor(x), method = cor_opt))

# Plot the resulting clustering tree (dendrogram)
pdf(file.path(diro, "03.tree.pdf"), width = 12, height = 9)
#sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE, hang = 0.04)
dev.off()

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
              deepSplit = 2, pamRespectsDendro = FALSE,
              minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file.path(diro, "04.tree.pdf"), width = 12, height = 9)
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                  dendroLabels = FALSE, hang = 0.03,
                  addGuide = TRUE, guideHang = 0.05,
                  main = "Gene dendrogram and module colors")
dev.off()
wg_clus1 = dynamicMods; wg_cols1 = dynamicColors

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf(file.path(diro, "05.module.cluster.pdf"), width = 7, height = 6)
#sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
dev.off()

#sizeGrWindow(12, 9)
pdf(file.path(diro, "06.pdf"), width = 12, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs
wg_clus2 = moduleLabels; wg_cols2 = mergedColors
wg_tree = geneTree

save(wg_tree, wg_clus1, wg_cols1, wg_clus2, wg_cols2, TOM, file = file.path(diro, "x.RData"))

         