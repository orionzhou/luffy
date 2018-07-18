require(plyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(RColorBrewer)
require(reshape2)
require(seqinr)
require(bios2mds)
source("Location.R")
source("comp.fun.R")

dirw = file.path(Sys.getenv("misc3"), "comp.ortho.hm")

## copy HM***_HM101/51_ortho/01.ortho.tbl to comp.ortho.hm/01_syn_ortho
for (qname in qnames_15) {
  f1 = sprintf("%s/%s_HM101/51_ortho/01.ortho.tbl", Sys.getenv("misc3"), qname)
  f2 = sprintf("%s/01_syn_ortho/%s.tbl", dirw, qname)
  system(sprintf("cp %s %s",f1, f2))
}

## create syn-ortho matrix
fg = file.path(tcfg$gdir, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:5,16)]

to = data.frame(idx = 1:nrow(tg), id = tg$id, stringsAsFactors = F)
colnames(to)[2] = tname
for (qname in qnames_15) {
  fgq = file.path(Sys.getenv("genome"), qname, "51.gtb")
  tgq = read.table(fgq, sep = "\t", header = T, as.is = T)
  qids = tgq$id

  fi = sprintf("%s/01_syn_ortho/%s.tbl", dirw, qname)
  ti = read.table(fi, sep = "\t", header = T, as.is = T)[,1:5]
  qidss = ti$qid[ti$qid != '']
  
  ti2 = ti[, c('tid','qid')]
	colnames(ti2) = c(tname, qname)

  to = merge(to, ti2, by = tname, all = T)
  cat(qname, sum(! qidss %in% qids), nrow(to), "\n")
}
to = to[order(to$idx), c(2,1,3:ncol(to))]

norgs = apply(to, 1, function(x) sum(x[-1] != ''))
table(norgs)

## run comp.ortho.score.R to create 05_score/xx.tbl
## compute pairwise aln similarity (HM101-based)
fi = file.path(dirw, "01.ids.tbl")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

fg = file.path(tcfg$gdir, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:5,16)]

tg = tg[order(tg$chr, tg$beg, tg$end),]
tg = cbind(tg, idx = 1:nrow(tg))
colnames(tg)[6] = 'fam'
tg = rename_genefam(tg)

to = tg[,1:2]
colnames(to)[1] = tname
for (qname in qnames_15) {
  fs = sprintf("%s/11_score/%s.tbl", dirw, qname)
  tt = read.table(fs, sep = "\t", header = T, as.is = T)
  tt = cbind(tt, cvg = 1 - (tt$qgap + tt$tgap) / tt$len, sim = tt$mat / (tt$mat + tt$mis))
  tt2 = tt[tt$cvg >= 0.5 & tt$sim >= 0.4, c('tid','sim')]
  colnames(tt2) = c(tname, qname)

  to = merge(to, tt2, by = tname, all = T)
  cat(qname, nrow(to), nrow(tt2), "\n")
}
to = to[,-2]
td = merge(ti[,1:2], to, by = tname)
stopifnot(nrow(td)==nrow(ti))
td = td[order(td$idx), c(2,1,3:ncol(td))]

norgs = apply(td, 1, function(x) sum(!is.na(x[-1])))
table(norgs)

fo = file.path(dirw, "12.score.tbl")
write.table(td, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

# filter syn-ortho matrix
sum(ti=='')
sum(is.na(td))
for (org in qnames_15) {
	idxs = which(is.na(td[,org]))
	ti[idxs, org] = ''
}
sum(ti=='')

fo = file.path(dirw, "05.ids.tbl")
write.table(ti, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

### add insertion to syn-ortho-matrix
fg = file.path(tcfg$gdir, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:5,16)]

fi = file.path(dirw, "05.ids.tbl")
ti = read.table(fi, sep = "\t", header = T, as.is = T)



### compute MPPD
fg = file.path(tcfg$gdir, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:5,16)]

fi = file.path(dirw, "05.ids.tbl")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

tos = ti[,c('idx', 'HM101', qnames_12)]
norgs = apply(tos, 1, function(x) sum(x[-1] != ''))
table(norgs)
toss = tos[norgs >= 10,]
nrow(toss)
fo = file.path(dirw, "21.ids.tbl")
write.table(toss, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')
# run comp.ortho.hm.pl create multiple-sequence alignment 

fi = file.path(dirw, "21.ids.tbl")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

tm = data.frame(id = tg$id, mpd = NA, stringsAsFactors = F)
for (idx in ti$idx) {
	fa = sprintf("%s/25_aln/%s.fas", dirw, idx)
	aln = read.alignment(fa, format = 'fasta')
	mat = dist.alignment(aln, matrix='similarity')
	tm$mpd[idx] = mean(mat)
	#aln = import.fasta(fa)
  #dis = mat.dis(aln, aln)
	#tm$mpd[idx] = mean(dis[lower.tri(dis)])
}
fo = file.path(dirw, "29.mpd.tbl")
write.table(tm, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

#### plot gene-fam syn-ortho view
fg = file.path(tcfg$gdir, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:5,16)]

tg = tg[order(tg$chr, tg$beg, tg$end),]
tg = cbind(tg, idx = 1:nrow(tg))
colnames(tg)[6] = 'fam'
tg = rename_genefam(tg)

fs = file.path(dirw, "12.score.tbl")
ts = read.table(fs, sep = "\t", header = T, as.is = T)
norgs = apply(ts, 1, function(x) sum(!is.na(x[-1])))
table(norgs)

ids_nomissing = ts$HM101[norgs >= 10]

## all pfam-misc
tgs = tg[tg$chr != 'chrU' & tg$id %in% ids_nomissing & !tg$fam %in% c("Unknown", "TE"),]
tgs = tgs[order(tgs$idx),]
tgs$idx = 1:nrow(tgs)

tx = ddply(tgs, .(chr), summarise, beg = min(idx), end = max(idx))
tx = tx[tx$end > tx$beg,]

tw = ts[ts$HM101 %in% tgs$id,]
colnames(tw)[1] = 'id'
tl = reshape(tw, direction = 'long', varying = list(2:16), idvar = 'id', timevar = 'org', v.names = 'score', times = colnames(tw)[2:16])

to = merge(tl, tgs[,c('id','idx')], by = 'id')
to$org = factor(to$org, levels = rev(qnames_15))

pb <- ggplot(to) +
  geom_tile(aes(x = idx, y = org, fill = 100*score), height = 0.8) +
  theme_bw() + 
  scale_x_continuous(name = '', limits = c(0, max(to$idx)+1), expand=c(0.01, 0), breaks = floor((tx$beg+tx$end)/2), labels = tx$chr) +
  scale_y_discrete(expand = c(0.02, 0), name = '') +
  scale_fill_gradient(name = 'sequence similarity with HM101 ortholog', space = "Lab", low = 'firebrick1', high = 'dodgerblue', na.value = 'grey50') +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0, 1), legend.title = element_text(size = 8), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill=NA, size=0)) +
  theme(plot.margin = unit(c(0,0.5,0,0), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid")) +
  annotate('segment', x = tx$beg, xend = tx$end, y = 0.3, yend = 0.3, size = 0.5, color = 'darkgreen') +
  annotate('segment', x = tx$beg, xend = tx$beg, y = 0.25, yend = 0.45, size = 0.5, color = 'darkgreen') +
  annotate('segment', x = tx$end, xend = tx$end, y = 0.25, yend = 0.45, size = 0.5, color = 'darkgreen')

fo = file.path(dirw, "13.score.pdf")
ggsave(pb, filename = fo, width = 8, height = 4)

## plot selected fams
famsp = c("Zinc-Finger", "CRP:NCR", "NBS-LRR")
pname = 'fams0'

plots = list()
for (i in 1:length(famsp)) {
  fam = famsp[i]

  tgs = tg[tg$chr != 'chrU' & tg$id %in% ids_nomissing & tg$fam == fam,]
  tgs = tgs[order(tgs$idx),]
  tgs$idx = 1:nrow(tgs)

  tx = ddply(tgs, .(chr), summarise, beg = min(idx), end = max(idx))
  tx = tx[tx$end > tx$beg,]

  tw = ts[ts$HM101 %in% tgs$id, -1]
  colnames(tw)[1] = 'id'
  tl = reshape(tw, direction = 'long', varying = list(2:16), idvar = 'id', timevar = 'org', v.names = 'score', times = colnames(tw)[2:16])

  to = merge(tl, tgs[,c('id','idx')], by = 'id')
  to$org = factor(to$org, levels = rev(qnames_15))

p1 <- ggplot(to) +
  geom_tile(aes(x = idx, y = org, fill = 100*score), height = 0.8) +
  theme_bw() + 
  scale_x_continuous(name = '', limits = c(0, max(to$idx)+1), expand=c(0.01, 0), breaks = floor((tx$beg+tx$end)/2), labels = tx$chr) +
  scale_y_discrete(expand = c(0.02, 0), name = '') +
  scale_fill_gradient(name = 'sequence similarity with HM101 ortholog', space = "Lab", low = 'firebrick1', high = 'dodgerblue', na.value = "white") +
  theme(plot.margin = unit(c(1,0.5,0,0), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid")) +
  annotate('segment', x = tx$beg, xend = tx$end, y = 0.3, yend = 0.3, size = 0.5, color = 'darkgreen') +
  annotate('segment', x = tx$beg, xend = tx$beg, y = 0.25, yend = 0.45, size = 0.5, color = 'darkgreen') +
  annotate('segment', x = tx$end, xend = tx$end, y = 0.25, yend = 0.45, size = 0.5, color = 'darkgreen')

  if(i == 1) {
    p1 <- p1 + theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0, 1), legend.title = element_text(size = 8), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "line"))
  } else {
    p1 <- p1 + theme(legend.position = 'none')
  }
  plots[[fam]] = p1
}

numrow = length(famsp)

fo = sprintf("%s/13.score.%s.pdf", dirw, pname)
pdf(file = fo, width = 6, height = 2.1*numrow+0.5, bg = 'transparent')
grid.newpage()
pushViewport(viewport(layout = grid.layout(numrow, 1, height = c(2.5,rep(2, numrow-1)))))

dco = data.frame(x = 1:numrow, y = rep(1,numrow), lab = paste(LETTERS[1:numrow], famsp, sep = ": "))
for (i in 1:nrow(dco)) {
  x = dco$x[i]; y = dco$y[i]; lab = dco$lab[i]
  print(plots[[famsp[i]]], vp = viewport(layout.pos.row = x, layout.pos.col = y))
  grid.text(lab, x = 0, y = unit(1, 'npc'), just = c('left', 'top'), 
    gp =  gpar(col = "black", fontface = 2, fontsize = 14),
    vp = viewport(layout.pos.row = x, layout.pos.col = y))
}
dev.off()


