require(ape)

#### plot NBS-LRR tree (with label)
dir = file.path(Sys.getenv('misc2'), 'nbs/mt_40')
fi = file.path(dir, "11.phb")
fo = file.path(dir, "11.png")

fm = file.path(dir, "../mt_10/13.flt.gal")
tm = read.table(fm, header = T, sep = "\t", as.is = T, quote = "")[, c(2,7)]

fc = file.path(dir, "../mt_10/clade.tbl")
tc = read.table(fc, header = T, sep = "\t", as.is = T, quote = "")[, c(1,4,6:7)]

tmc = merge(tm, tc, by.x = 'qId', by.y = 'tid')

tree = read.tree(fi)

pick_clade <- function(id, df) {
  idx = which(df$tId == id)
  ifelse(length(idx) == 0, '', df$clade[idx[1]])
}
labels = tree$tip.label
clades = sapply(labels, pick_clade, tmc)
tree$tip.label = paste(labels, clades, sep = " ")

  scores = as.numeric(tree$node.label)
  if(mean(scores, na.rm = TRUE) > 1) { scores = scores / 1000 }
  node.labels.bg = rep('white', tree$Nnode)
  node.labels.bg[scores >= 0.95] = 'black'
  node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

png(filename = fo, width = 800, height = 7000, units = 'px')
plot(tree, show.node.label = F, show.tip.label = T, font = 1, 
  label.offset = 0.01, no.margin = T, cex = 0.78)
nodelabels(pch = 22, bg = node.labels.bg)
add.scale.bar(lcol = 'black')
dev.off()

tt = data.frame(id = labels, clade = clades)
write.table(tt, file.path(dir, '15.clade.tbl'), row.names = F, col.names = T,
  sep = "\t", quote = F)

#### plot NBS-LRR tree (simple)
dir = file.path(Sys.getenv('misc2'), 'nbs/mt_40')
fi = file.path(dir, "26.phb")
fo = file.path(dir, "26.png")
tree = read.tree(fi)

png(filename = fo, width = 800, height = 7000, units = 'px')
plot(tree, show.tip.label = T, font = 1, 
  label.offset = 0.01, no.margin = T, cex = 0.78)
add.scale.bar(lcol = 'black')
dev.off()