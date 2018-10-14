require(ape)
require(seqinr)
require(RColorBrewer)

dirw = file.path(Sys.getenv('misc2'), 'bact')

### extract sequences at the class level
fc = file.path(dirw, "classes.txt")
tc = read.table(fc, header = F, sep = "\t", as.is = T, quote= "")
cls = unique(tc$V1)

fl = file.path(dirw, "11.fas.tsv")
tl = read.table(fl, header = F, sep = "\t", as.is = T, quote = "")
colnames(tl) = c("id", "taxa")
labs = sapply(strsplit(tl$taxa, split=";"), "[", 6)
idxs = which(labs %in% cls)
stopifnot(length(cls) == length(idxs))

fsi = file.path(dirw, "11.fas")
seqs <- read.fasta(fsi, as.string=T)
fso = file.path(dirw, "12.fas")
seqss = seqs[idxs]
write.fasta(sequences = seqss, names = names(seqss), nbchar = 80, file.out = fso)
# phy.py 12.fas 23

### plot tree
ft = file.path(dirw, "23.nwk")
tree = read.tree(ft)

fl = file.path(dirw, "11.fas.tsv")
tl = read.table(fl, header = F, sep = "\t", as.is = T, quote = "")
colnames(tl) = c("id", "taxa")
genus = sapply(strsplit(tl$taxa, split=";"), "[", 4)
#genus = sapply(strsplit(tl$taxa, split=";"), "[", 6)
#genus = sapply(strsplit(genus, split=" "), "[", 1)
tl = cbind.data.frame(tl, genus = as.character(genus), stringsAsFactors = F)

ids = tree$tip.label
ntip = length(ids)
tl = tl[tl$id %in% ids,]
stopifnot(nrow(tl) == length(ids))
tl = tl[match(ids, tl$id),]

rootlab = 'Opitutales' #Alterococcus

tree$tip.label = tl$genus
tree = root(tree, which(tree$tip.label==rootlab)) 
#tree = root(tree, node=tree$edge[which(tree$edge[,2]==which.edge(tree,rootlab)),1])
tree = rotate(tree, 23)

#  genus1 = c('Bifidobacterium', 'Microbispora', 'Actinomadura', 'Kocuria', 'Clavibacter', 'Planomonospora')
#  genus2 = c('Actinoplanes')
#  genus3 = c('Streptomyces')
  
  tip.cols = rep('black', length(ids))
  tip.cols[tl$genus == rootlab] = 'seagreen'
#  tip.cols[tl$genus %in% genus1] = brewer.pal(9, "Set1")[3]
#  tip.cols[tl$genus %in% genus2] = brewer.pal(9, "Set1")[2]
#  tip.cols[tl$genus %in% genus3] = brewer.pal(9, "Set1")[1]
	
	labs = rep("", length(ids))
	labs[tl$genus == 'Streptosporangiales'] = 3
	labs[tl$genus == 'Micrococcales'] = 2
	labs[tl$genus == 'Micromonosporales'] = 3
	labs[tl$genus == 'Streptomycetales'] = 10
	labs[tl$genus == 'Bifidobacteriales'] = 1
  
  scores = as.numeric(tree$node.label)
  node.bg = rep('white', tree$Nnode)
  node.bg[scores >= 0.9] = 'black'
  node.bg[scores >= 0.8 & scores < 0.9] = 'grey'
  tree$node.label = sprintf("%.02f", as.numeric(tree$node.label))
  
  fo = file.path(dirw, "25.pdf")
  pdf(fo, width = 5, height = 6)
  plot(tree, type = 'phylogram', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols, label.offset = 0.01, 
    no.margin = T, cex = 0.8, x.lim = 1.1, y.lim = ntip+2, align.tip.label = T)
  #tiplabels(labs, frame = 'none', adj = 1)
  nodelabels(pch = 22, bg = node.bg, cex = 0.5)
  text(1, 1:ntip, labels = labs, cex = 0.9)
  text(1, ntip+1.2, labels = '# of Lantibiotics', cex = 0.8)
  text(1, ntip+0.7, labels = 'Identified', cex = 0.8)
  add.scale.bar(x = 0, y = ntip+1, lcol = 'black')
  
#  legend(0.58, 29, title = "# of Lantibiotics Characterized", legend = c("1 ", "3", "9"), fill = brewer.pal(9, "Set1")[3:1], xjust = 0.5, cex=0.8)
  dev.off()

##### plot tree - family level
ft = file.path(dirw, "56.nwk")
tree = read.tree(ft)

fl = file.path(dirw, "51.fas.tsv")
tl = read.table(fl, header = F, sep = "\t", as.is = T, quote = "")
colnames(tl) = c("id", "taxa")
genus = sapply(strsplit(tl$taxa, split=";"), "[", 5)
tl = cbind.data.frame(tl, genus = as.character(genus), stringsAsFactors = F)

ids = tree$tip.label
ntip = length(ids)
tl = tl[tl$id %in% ids,]
stopifnot(nrow(tl) == length(ids))
tl = tl[match(ids, tl$id),]

rootlab = 'Opitutaceae' #Alterococcus

tree$tip.label = tl$genus
tree = root(tree, which(tree$tip.label==rootlab)) 
#tree = root(tree, node=tree$edge[which(tree$edge[,2]==which.edge(tree,rootlab)),1])
tree = rotate(tree, 74)
  
  tip.cols = rep('black', length(ids))
  tip.cols[tl$genus == rootlab] = 'dodgerblue'
  plabs = c('Streptosporangiaceae', 'Micrococcaceae', 'Micromonosporaceae', 'Streptomycetaceae', 'Bifidobacteriaceae')
  tip.cols[tl$genus %in% plabs] = 'brown1'
	
	labs = rep("", length(ids))
	labs[tl$genus == 'Streptosporangiaceae'] = 3
	labs[tl$genus == 'Micrococcaceae'] = 2
	labs[tl$genus == 'Micromonosporaceae'] = 3
	labs[tl$genus == 'Streptomycetaceae'] = 10
	labs[tl$genus == 'Bifidobacteriaceae'] = 1
  
  scores = as.numeric(tree$node.label)
  node.bg = rep('white', tree$Nnode)
  node.bg[scores >= 0.9] = 'black'
  node.bg[scores >= 0.8 & scores < 0.9] = 'grey'
  tree$node.label = sprintf("%.02f", as.numeric(tree$node.label))
  
  fo = file.path(dirw, "59.pdf")
  pdf(fo, width = 8, height = 10)
  plot(tree, type = 'phylogram', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols, label.offset = 0.01, 
    no.margin = T, cex = 0.8, x.lim = 1.5, y.lim = ntip+2, align.tip.label = T)
  #tiplabels(labs, frame = 'none', adj = 1)
  nodelabels(pch = 22, bg = node.bg, cex = 0.5)
  text(1.4, 1:ntip, labels = labs, cex = 0.9)
  text(1.4, ntip+2.3, labels = '# of Lantibiotics', cex = 0.8)
  text(1.4, ntip+1.2, labels = 'Identified', cex = 0.8)
  add.scale.bar(x = 0, y = ntip+1, lcol = 'black')
  dev.off()

### plot tree
ft = file.path(dirw, "66.nwk")
tree = read.tree(ft)

ids = tree$tip.label
ntip = length(ids)

rootlab = 'Opitutaceae' #Alterococcus

tree = root(tree, which(tree$tip.label==rootlab)) 
#tree = root(tree, node=tree$edge[which(tree$edge[,2]==which.edge(tree,rootlab)),1])
tree = rotate(tree, 56)
  
  tip.cols = rep('black', length(ids))
  tip.cols[ids == rootlab] = 'dodgerblue'
  plabs = c('Streptosporangiaceae', 'Micrococcaceae', 'Micromonosporaceae', 'Streptomycetaceae', 'Bifidobacteriaceae')
#  tip.cols[tl$genus %in% plabs] = 'brown1'
	
	labs = rep("", length(ids))
	labs[ids == 'Streptosporangiaceae'] = 3
	labs[ids == 'Micrococcaceae'] = 2
	labs[ids == 'Micromonosporaceae'] = 3
	labs[ids == 'Streptomycetaceae'] = 10
	labs[ids == 'Bifidobacteriaceae'] = 1
  
  scores = as.numeric(tree$node.label)
  node.bg = rep('white', tree$Nnode)
  node.bg[scores >= 0.8] = 'black'
  #node.bg[scores >= 0.8 & scores < 0.9] = 'grey'
  tree$node.label = sprintf("%.02f", as.numeric(tree$node.label))
  
  fo = file.path(dirw, "69.pdf")
  pdf(fo, width = 6, height = 8)
  plot(tree, type = 'phylogram', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols, label.offset = 0.01, , font = 1,
    no.margin = T, cex = 0.75, x.lim = 1.9, y.lim = ntip+2, align.tip.label = T)
  #tiplabels(labs, frame = 'none', adj = 1)
  nodelabels(pch = 22, bg = node.bg, cex = 0.7)
  text(1.8, 1:ntip, labels = labs, cex = 0.8)
  text(1.8, ntip+2.3, labels = '# of Lantibiotics', cex = 0.75)
  text(1.8, ntip+1.2, labels = 'Identified', cex = 0.75)
  add.scale.bar(x = 0, y = ntip+2, lcol = 'black')
  dev.off()
  