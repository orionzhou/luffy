require(ape)

f_ann = "/home/youngn/zhoux379/data/misc3/hapmap/31_phylogeny/mt_label.tbl"
ann = read.table(f_ann, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

plot_mt_tree_2 <- function(fi, fo, ann, opt) {
  tree = read.tree(fi)

  label_id = tree$tip.label
  idx_HM101 = which(label_id == 'HM101')
  idx_acc26 = which(label_id %in% get_mt_ids('acc26'))

  tip.color = rep('black', length(tree$tip.label))
  tip.color[idx_acc26] = 'royalblue'
  tip.color[idx_HM101] = 'red'
  
  df1 = data.frame(id=tree$tip.label)
  df2 = merge(df1, ann, by="id", all.x=TRUE)
  label_origin = paste(df2$country, df2$category, sep=" | ")
  tree$tip.label = label_origin

  scores = tree$node.label
  node.labels.bg = rep('white', tree$Nnode)
  node.labels.bg[scores >= 0.95] = 'black'
  node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

  if(opt == 'acc26') {
    cex1 = 1
    cex2 = 1.2
    label.offset = 0.07
    xmax = 0.5
  } else if (opt == 'acc31') {
    cex1 = 1
    cex2 = 1.2
    label.offset = 0.07
    xmax = 0.75
  } else {
    stop(cat("unknown opt: ", opt, "\n", sep=""))
  }

  png(filename=fo, width=600, height=600, units='px')
  plot(tree, show.node.label=FALSE, show.tip.label=TRUE, tip.color=tip.color, font=3, no.margin=TRUE, 
    cex=cex1, label.offset=label.offset, x.lim=xmax)
  tiplabels(label_id, adj=0, bg=NA, font=2, frame='none', col=tip.color, cex=cex2)
  nodelabels(pch = 21, bg=node.labels.bg)
  add.scale.bar()
  dev.off()
}

### plot acc31 tree
opt = 'acc31'
dir = file.path(DIR_Data, "misc3/hapmap_mt35/31_phylogeny", opt)
for(i in 1:8) {
  fi = sprintf("%s/22_phyml/chr%d.nwk", dir, i)
  fo = sprintf("%s/22_phyml/chr%d.png", dir, i)
  plot_mt_tree(fi, fo, ann)
}
fi = file.path(dir, "19/04.nwk")
fo = file.path(dir, "19/04.png")
plot_mt_tree(fi, fo)

### plot acc56 tree
opt = 'acc56'
dir = file.path(DIR_Data, "misc3/hapmap_mt35/31_phylogeny", opt)
for(i in 1:8) {
  fi = sprintf("%s/22_phyml/chr%d.nwk", dir, i)
  fo = sprintf("%s/22_phyml/chr%d.png", dir, i)
  plot_mt_tree(fi, fo, ann)
}

### plot deepseq tree
opt = 'deepseq'
reg = "chr5"
dir = file.path(DIR_Data, "misc3/hapmap/31_phylogeny", opt)
fi = sprintf("%s/22_phyml/%s.nwk", dir, reg)
fo = sprintf("%s/22_phyml/%s.png", dir, reg)
#fi = sprintf("%s/21_phynj/%s.phb", dir, reg)
#fo = sprintf("%s/21_phynj/%s.png", dir, reg)
  tree = read.tree(fi)

  group1 = c("HM101", "HM056", "HM058", "HM117", "HM125")
  group2 = c("HM340", "HM324", "HM018", "HM022-I", "HM017-I")
  group3 = c("HM034", "HM129", "HM060", "HM095", "HM185")
  grouph = c("HM034", "HM056", "HM340")
  
  labels = tree$tip.label
  tip.color = rep('black', length(tree$tip.label))
  tip.color[which(labels %in% group1)] = 'red'
  tip.color[which(labels %in% group2)] = 'forestgreen'
  tip.color[which(labels %in% group3)] = 'dodgerblue'
  
  font = rep(1, length(tree$tip.label))
  font[which(labels %in% grouph)] = 2
  font[which(labels == "HM101")] = 4

  df1 = data.frame(idx = 1 : length(labels), id = labels)
  df2 = merge(df1, ann, by = "id", all.x = T)
  df3 = df2[order(df2$idx), ]
  labelsn = as.character(df3$id)
  for (i in 1:nrow(df3)) {
    id = df3$id[i]
    country = df3$country[i]
    if(!is.na(country) & country != "") { 
      labelsn[i] = paste(df3$id[i], df3$country[i], sep = " | ")
    }
  }
  tree$tip.label = labelsn

  scores = as.numeric(tree$node.label)
  if(mean(scores, na.rm = TRUE) > 1) { scores = scores / 1000 }
  node.labels.bg = rep('white', tree$Nnode)
  node.labels.bg[scores >= 0.95] = 'black'
  node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

  png(filename=fo, width=600, height=600, units='px')
  plot(tree, show.node.label = F, show.tip.label = T, font = font, 
    tip.color = tip.color, label.offset = 0.005, no.margin = T, cex = 1.2)
  nodelabels(pch = 22, bg = node.labels.bg)
  add.scale.bar(x = 0.02, y = 20, lcol = 'black')
  dev.off()

### plot ingroup tree
opt = 'ingroup'
reg = "chr5"
dir = file.path(DIR_Data, "misc3/hapmap/31_phylogeny", opt)
fi = sprintf("%s/22_phyml/%s.nwk", dir, reg)
fo = sprintf("%s/22_phyml/%s.png", dir, reg)
#fi = sprintf("%s/21_phynj/%s.phb", dir, reg)
#fo = sprintf("%s/21_phynj/%s.png", dir, reg)
  tree = read.tree(fi)

group1 = c(
  "HM058", "HM125", "HM056", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022-I", "HM324", "HM340"
)
group2 = c("HM101", "HM034", "HM056", "HM340")
grouph = c("HM101")
  
  labels = tree$tip.label
  tip.color = rep('black', length(tree$tip.label))
  tip.color[which(labels %in% group1)] = 'red'
  tip.color[which(labels %in% group2)] = 'forestgreen'
#  tip.color[which(labels %in% group3)] = 'dodgerblue'
#  tip.color[which(labels %in% group4)] = 'purple'
  
  font = rep(1, length(tree$tip.label))
  font[which(labels %in% grouph)] = 2
  font[which(labels == "HM101")] = 4

  df1 = data.frame(idx = 1 : length(labels), id = labels)
  df2 = merge(df1, ann, by = "id", all.x = T)
  df3 = df2[order(df2$idx), ]
  labelsn = as.character(df3$id)
  for (i in 1:nrow(df3)) {
    id = df3$id[i]
    country = df3$country[i]
    if(!is.na(country) & country != "") { 
      labelsn[i] = paste(df3$id[i], df3$country[i], sep = " | ")
    }
  }
  tree$tip.label = labelsn

  scores = as.numeric(tree$node.label)
  if(mean(scores, na.rm = TRUE) > 1) { scores = scores / 1000 }
  node.labels.bg = rep('white', tree$Nnode)
  node.labels.bg[scores >= 0.95] = 'black'
  node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

  png(filename=fo, width=1000, height=4000, units='px')
  plot(tree, show.node.label = F, show.tip.label = T, font = font, 
    tip.color = tip.color, label.offset = 0.005, no.margin = T, cex = 1.2)
  nodelabels(pch = 22, bg = node.labels.bg)
  add.scale.bar(x = 0.02, y = 20, lcol = 'black')
  dev.off()

### plot RIL1 tree
opt = 'ril1'
reg = "chr5"
dir = file.path(DIR_Data, "misc3/hapmap/31_phylogeny", opt)
fi = sprintf("%s/22_phyml/%s.nwk", dir, reg)
fo = sprintf("%s/22_phyml/%s.png", dir, reg)

  tree = read.tree(fi)

  group1 = c("HM004", "HM005", "HM006", "HM017-I", "HM018", "HM019", "HM026", "HM027", "HM028")
  group2 = c("HM101")
  
  labels = tree$tip.label
  tip.color = rep('black', length(tree$tip.label))
  tip.color[which(labels %in% group1)] = 'dodgerblue'
  tip.color[which(labels %in% group2)] = 'red'
#  tip.color[which(labels %in% group3)] = 'forestgreen'
  
  font = rep(2, length(tree$tip.label))
#  font[which(labels %in% grouph)] = 2

  df1 = data.frame(idx = 1 : length(labels), id = labels)
  df2 = merge(df1, ann, by = "id", all.x = T)
  df3 = df2[order(df2$idx), ]
  labelsn = as.character(df3$id)
  for (i in 1:nrow(df3)) {
    id = df3$id[i]
    country = df3$country[i]
    if(!is.na(country) & country != "") { 
      labelsn[i] = paste(df3$id[i], df3$country[i], sep = " | ")
    }
  }
  tree$tip.label = labelsn

  scores = as.numeric(tree$node.label)
  if(mean(scores, na.rm = TRUE) > 1) { scores = scores / 1000 }
  node.labels.bg = rep('white', tree$Nnode)
  node.labels.bg[scores >= 0.95] = 'black'
  node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

  png(filename=fo, width=600, height=600, units='px')
  plot(tree, show.node.label = F, show.tip.label = T, font = font, x.lim = 0.6,
    tip.color = tip.color, label.offset = 0.005, no.margin = T, cex = 1.3)
  nodelabels(pch = 22, bg = node.labels.bg)
  add.scale.bar(x = 0.02, y = 13, lcol = 'black')
  dev.off()

### plot RIL2 tree
opt = 'ril2'
reg = "chr5"
dir = file.path(DIR_Data, "misc3/hapmap/31_phylogeny", opt)
fi = sprintf("%s/22_phyml/%s.nwk", dir, reg)
fo = sprintf("%s/22_phyml/%s.png", dir, reg)

  tree = read.tree(fi)

  group1 = c("HM004", "HM005", "HM006", "HM017-I", "HM018", "HM019", "HM026", "HM027", "HM028")
  group2 = c("HM101")
  
  labels = tree$tip.label
  tip.color = rep('black', length(tree$tip.label))
  tip.color[which(labels %in% group1)] = 'dodgerblue'
  tip.color[which(labels %in% group2)] = 'red'
#  tip.color[which(labels %in% group3)] = 'forestgreen'
  
  font = rep(2, length(tree$tip.label))
#  font[which(labels %in% grouph)] = 2

  df1 = data.frame(idx = 1 : length(labels), id = labels)
  df2 = merge(df1, ann, by = "id", all.x = T)
  df3 = df2[order(df2$idx), ]
  labelsn = as.character(df3$id)
  for (i in 1:nrow(df3)) {
    id = df3$id[i]
    country = df3$country[i]
    if(!is.na(country) & country != "") { 
      labelsn[i] = paste(df3$id[i], df3$country[i], sep = " | ")
    }
  }
  tree$tip.label = labelsn

  scores = as.numeric(tree$node.label)
  if(mean(scores, na.rm = TRUE) > 1) { scores = scores / 1000 }
  node.labels.bg = rep('white', tree$Nnode)
  node.labels.bg[scores >= 0.95] = 'black'
  node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

  png(filename = fo, width = 600, height = 800, units = 'px')
  plot(tree, show.node.label = F, show.tip.label = T, font = font, x.lim = 0.55,
    tip.color = tip.color, label.offset = 0.005, no.margin = T, cex = 1.2)
  nodelabels(pch = 22, bg = node.labels.bg)
  add.scale.bar(x = 0.02, y = 22, lcol = 'black')
  dev.off()

### plot pan16-denovo tree
dir = file.path(DIR_Data, "misc3/phy_denovo")
fi = file.path(dir, "12.nwk")
fo = file.path(dir, "12.png")
tree = read.tree(fi)

group1 = c("HM101", "HM056", "HM058", "HM117", "HM125")
group2 = c("HM340", "HM324", "HM018", "HM022", "HM017")
group3 = c("HM034", "HM129", "HM060", "HM095", "HM185")
grouph = c("HM034", "HM056", "HM340")

labels = tree$tip.label
tip.color = rep('black', length(tree$tip.label))
tip.color[which(labels %in% group1)] = 'red'
tip.color[which(labels %in% group2)] = 'forestgreen'
tip.color[which(labels %in% group3)] = 'dodgerblue'

font = rep(1, length(tree$tip.label))
font[which(labels %in% grouph)] = 2
font[which(labels == "HM101")] = 4

df1 = data.frame(idx = 1 : length(labels), id = labels)
df2 = merge(df1, ann, by = "id", all.x = T)
df3 = df2[order(df2$idx), ]
labelsn = as.character(df3$id)
for (i in 1:nrow(df3)) {
  id = df3$id[i]
  country = df3$country[i]
  if(!is.na(country) & country != "") { 
    labelsn[i] = paste(df3$id[i], df3$country[i], sep = " | ")
  }
}
tree$tip.label = labelsn

scores = as.numeric(tree$node.label)
if(mean(scores, na.rm = TRUE) > 1) { scores = scores / 1000 }
node.labels.bg = rep('white', tree$Nnode)
node.labels.bg[scores >= 0.95] = 'black'
node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

png(filename=fo, width=600, height=600, units='px')
plot(tree, show.node.label = F, show.tip.label = T, font = font, 
  tip.color = tip.color, label.offset = 0.005, no.margin = T, cex = 1.2)
nodelabels(pch = 22, bg = node.labels.bg)
add.scale.bar(x = 0.02, y = 10, lcol = 'black')
dev.off()

### plot ingroup + HM018BC + HM022B tree
dirw = file.path(DIR_Data, "misc3/phy_finderror")
fi = file.path(dirw, "52.nwk")
fo = file.path(dirw, "52.pdf")
tree = read.tree(fi)

group1 = c("HM018", "HM018B", "HM018C")
group2 = c("HM022-I", "HM022B")
group3 = c("HM017-I", "HM017B")
grouph = c("HM101", "HM340")

labels = tree$tip.label
tip.color = rep('black', length(tree$tip.label))
tip.color[which(labels %in% group1)] = 'red'
tip.color[which(labels %in% group2)] = 'forestgreen'
tip.color[which(labels %in% group3)] = 'purple'
tip.color[which(labels %in% grouph)] = 'dodgerblue'

font = rep(1, length(tree$tip.label))
font[which(labels %in% grouph)] = 2
font[which(labels == "HM101")] = 4

df1 = data.frame(idx = 1 : length(labels), id = labels)
df2 = merge(df1, ann, by = "id", all.x = T)
df3 = df2[order(df2$idx), ]
labelsn = as.character(df3$id)
for (i in 1:nrow(df3)) {
  id = df3$id[i]
  country = df3$country[i]
  if(!is.na(country) & country != "") { 
    labelsn[i] = paste(df3$id[i], df3$country[i], sep = " | ")
  }
}
tree$tip.label = labelsn

scores = as.numeric(tree$node.label)
if(mean(scores, na.rm = TRUE) > 1) { scores = scores / 1000 }
node.labels.bg = rep('white', tree$Nnode)
node.labels.bg[scores >= 0.95] = 'black'
node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

pdf(fo, width = 8, height = 32)
plot(tree, show.node.label = F, show.tip.label = T, font = font, 
  tip.color = tip.color, label.offset = 0.005, no.margin = T, cex = 0.75)
nodelabels(pch = 22, bg = node.labels.bg)
add.scale.bar(x = 0.02, y = 50, lcol = 'black')
dev.off()

### plot pan16 tree
dir = file.path(DIR_Data, "misc1/phy.mt/pan16")
fi = file.path(dir, "31.nwk")
fo = file.path(dir, "32.pdf")
tree = read.tree(fi)

grouph = c("HM101", "HM340")
group2 = c("HM101", "HM034", "HM056", "HM340")

labs = tree$tip.label

font = rep(1, length(labs))
#font[which(labs %in% grouph)] = 4
tip.color = rep('black', length(labs))
tip.color[which(labs %in% grouph)] = 'dodgerblue'

label2.bg = rep('white', length(labs))
label2.bg[which(labs %in% group2)] = 'red'

df1 = data.frame(idx = 1:length(labs), id = labs)
df2 = merge(df1, ann, by = "id", all.x = T)
df3 = df2[order(df2$idx), ]
notes = as.character(df3$country)
notes[is.na(notes)] = ""
notes[labs == "HM101"] = paste("(A17)", notes[labs == "HM101"], by = " ")
notes[labs == "HM340"] = paste("(R108)", notes[labs == "HM340"], by = " ")
notes[labs == 'HM023'] = "Tunisia"
#notes = sprintf(paste("%-", max(nchar(notes)), "s", sep = ''), notes)

scores = as.numeric(tree$node.label)
if(mean(scores, na.rm = TRUE) > 1) { scores = scores / 1000 }
node.labels.bg = rep('white', tree$Nnode)
node.labels.bg[scores >= 0.95] = 'black'
node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

tree$tip.label = paste(labs, notes, sep = "       ")
#tree = root(tree, 4)
pdf(file = fo, width = 5, height = 6, bg = 'transparent')
plot(tree, show.node.label = F, show.tip.label = T, font = font, x.lim = 0.8,
  tip.color = tip.color, label.offset = 0.005, no.margin = T, cex = 0.9)
nodelabels(pch = 22, bg = node.labels.bg)
tiplabels(pch = 22, col = NA, bg = label2.bg, adj = 0.61, cex = 1.5)
#par(family = "Courier New")
add.scale.bar(x = 0, y = 10, lcol = 'black')
rect(0.08, 11.87, 0.10, 12.13, col = 'red', border = NA)
text(0.115, 12, labels = "RNA-Seq", cex = 0.75, adj = c(0, 0.5))
dev.off()

### plot pan16-expanded tree
dir = file.path(DIR_Data, "misc1/phy.mt/pan16x")
fi = file.path(dir, "31.nwk")
fo = file.path(dir, "32.pdf")
tree = read.tree(fi)

grouph = c("HM101")
group1 = c(
  "HM058", "HM125", "HM056", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340"
)
group2 = c("HM101", "HM034", "HM056", "HM340")

labs = tree$tip.label

font = rep(1, length(labs))
font[which(labs %in% grouph)] = 4
tip.color = rep('black', length(labs))
tip.color[which(labs %in% grouph)] = 'dodgerblue'

label1.bg = rep('white', length(labs))
label1.bg[which(labs %in% group1)] = 'red'

label2.bg = rep('white', length(labs))
label2.bg[which(labs %in% group2)] = 'forestgreen'

df1 = data.frame(idx = 1:length(labs), id = labs)
df2 = merge(df1, ann, by = "id", all.x = T)
df3 = df2[order(df2$idx), ]
notes = as.character(df3$country)
notes[is.na(notes)] = ""
#notes = sprintf(paste("%-", max(nchar(notes)), "s", sep = ''), notes)

scores = as.numeric(tree$node.label)
if(mean(scores, na.rm = TRUE) > 1) { scores = scores / 1000 }
node.labels.bg = rep('white', tree$Nnode)
node.labels.bg[scores >= 0.95] = 'black'
node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

tree$tip.label = paste(labs, notes, sep = "       ")
#tree = root(tree, 4)
pdf(file = fo, width = 6, height = 7.5, bg = 'transparent')
plot(tree, show.node.label = F, show.tip.label = T, font = font, x.lim = 0.67,
  tip.color = tip.color, label.offset = 0.005, no.margin = T, cex = 0.9)
nodelabels(pch = 22, bg = node.labels.bg)
tiplabels(pch = 22, col = NA, bg = label1.bg, adj = 0.573, cex = 1.5)
tiplabels(pch = 22, col = NA, bg = label2.bg, adj = 0.587, cex = 1.5)
#par(family = "Courier New")
add.scale.bar(x = 0, y = 16, lcol = 'black')
rect(0.1, 25, 0.113, 25.6, col = 'red', border = NA)
rect(0.1, 24, 0.113, 24.6, col = 'forestgreen', border = NA)
text(0.12, 25.3, labels = "ALLPATHS-LG", cex = 0.7, adj = c(0, 0.5))
text(0.12, 24.3, labels = "RNA-Seq", cex = 0.7, adj = c(0, 0.5))
dev.off()

### hapmap-denovo
dirw = file.path(DIR_Data, "misc1/phy.mt/hapmapdenovo")
fi = file.path(dirw, "31.nwk")
fo = file.path(dirw, "32.1.pdf")
tree = read.tree(fi)

grouph = c("HM101")
groupo = c("HM340")

labs = tree$tip.label

font = rep(1, length(labs))
font[which(labs %in% c(grouph, groupo))] = 4
tip.color = rep('black', length(labs))
tip.color[which(labs %in% grouph)] = 'dodgerblue'
tip.color[which(labs %in% groupo)] = 'forestgreen'

df1 = data.frame(idx = 1:length(labs), id = labs)
df2 = merge(df1, ann, by = "id", all.x = T)
df3 = df2[order(df2$idx), ]
notes = as.character(df3$country)
notes[is.na(notes)] = ""
#notes = sprintf(paste("%-", max(nchar(notes)), "s", sep = ''), notes)

scores = as.numeric(tree$node.label)
if(mean(scores, na.rm = TRUE) > 1) { scores = scores / 1000 }
node.labels.bg = rep('white', tree$Nnode)
node.labels.bg[scores >= 0.95] = 'black'
node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

tree$tip.label = paste(labs, notes, sep = "   ")
tree$tip.label = rep('', length(labs))
#tree = root(tree, 4)
pdf(file = fo, width = 5, height = 6, bg = 'transparent')
plot(tree, show.node.label = F, show.tip.label = T, font = font, x.lim = 0.6,
  tip.color = tip.color, label.offset = 0.005, no.margin = T, cex = 0.9)
nodelabels(pch = 22, bg = node.labels.bg)
#par(family = "Courier New")
add.scale.bar(x = 0, y = 16, lcol = 'black')
dev.off()

### compare sv phylogeny with chr5 phylogeny
  fi_chr5 = file.path(DIR_Data, "repo/mt_35/31_phylogeny", "acc26", "21_phynj/chr5.phb")
  fo_chr5 = file.path(DIR_Data, "repo/mt_35/31_phylogeny", "acc26", "21_phynj/chr5.png")
  tree = read.tree(fi_chr5)

  tip.text = tree$tip.label
  tip.bg.chr5 = rainbow(length(tree$tip.label))
  names(tip.bg.chr5) = tip.text

  scores = as.numeric(tree$node.label)
  if(mean(scores, na.rm=TRUE) > 1) { scores = scores / 1000 }
  node.bg = rep('white', tree$Nnode)
  node.bg[scores >= 0.95] = 'black'
  node.bg[scores >= 0.8 & scores < 0.95] = 'gray'

  png(filename=fo_chr5, width=500, height=500, units='px')
  plot(tree, show.node.label=FALSE, show.tip.label=FALSE, font=4, no.margin=TRUE, cex=1.2)
  tiplabels(text=tip.text, bg=tip.bg.chr5)
  nodelabels(pch=22, bg=node.bg)
  add.scale.bar()
  dev.off()

  fi_sv = file.path(DIR_Repo, "mt_35/40_sv/41_shared/73_phylogeny", "04.phb")
  fo_sv = file.path(DIR_Repo, "mt_35/40_sv/41_shared/73_phylogeny", "04.png")
  tree = read.tree(fi_sv)
  tree = rotate(tree, c("HM015", "HM003"))
#  tree = rotate(tree, c("HM101", "HM008"))
#  tree = rotate(tree, c("HM101", "HM007"))
  tree = rotate(tree, c("HM020", "HM021"))
  tree = rotate(tree, c("HM021", "HM028"))

  tip.text = tree$tip.label
  tip.bg.sv = c(tip.bg.chr5[tip.text])

  scores = as.numeric(tree$node.label)
  if(mean(scores, na.rm=TRUE) > 1) { scores = scores / 1000 }
  node.bg = rep('white', tree$Nnode)
  node.bg[scores >= 0.95] = 'black'
  node.bg[scores >= 0.8 & scores < 0.95] = 'gray'

  png(filename=fo_sv, width=500, height=500, units='px')
  plot(tree, show.node.label=FALSE, show.tip.label=FALSE, font=4, no.margin=TRUE, cex=1.2)
  tiplabels(text=tip.text, bg=tip.bg.sv)
  nodelabels(pch=22, bg=node.bg)
  add.scale.bar()
  dev.off()

  

