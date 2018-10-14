require(seqinr)
require(geiger)
require(igraph)
require(grid)
require(plyr)
require(ggplot2)
source("Align.R")
source(file.path(Sys.getenv('code'), 'r', 'clustertree.R'))


org = "HM340.AC"
dirw = file.path(Sys.getenv("misc2"), "bd.crp", org)
dir.create(dirw)
setwd(dirw)

f_gen = file.path(Sys.getenv("genome"), org, "51.gtb")
tg = read.table(f_gen, header = T, sep = "\t", as.is = T)[,c(1,3:5,15:17)]
tg = tg[order(tg$chr, tg$beg, tg$end),]
tg = cbind(idx = 1:nrow(tg), tg)

tg = tg[tg$cat2 %in% c("CRP0000-1030", "NCR", "CRP1600-6250"),]

tb = table(tg$cat3)
fams = names(tb[tb>=5])

tg2 = tg[tg$cat3 %in% fams,]

### generate protein / cds sequnces
#gtb2fas.pl -i 43.crp.gtb -o 43.crp.cds.fas -d 11_genome.fas -p cds
#gtb2fas.pl -i 43.crp.gtb -o 43.crp.pro.fas -d 11_genome.fas -p protein

### make protein / cds alignments
cdsseq <- read.fasta(file.path(Sys.getenv("genome"), org, "43.crp.cds.fas"), seqtype = "DNA", as.string = T, set.attributes = F)
proseq <- read.fasta(file.path(Sys.getenv("genome"), org, "43.crp.pro.fas"), seqtype = "AA", as.string = T, set.attributes = F)

dir.create("01.cds.seq"); dir.create("01.pro.seq"); dir.create("03.pro.aln")
for (fam in fams) {
  to = tg2[tg2$cat3 == fam,]
  ids = to$id
  
  cseq = cdsseq[ids]
  foc = sprintf("%s/01.cds.seq/%s.fas", dirw, fam)
  write.fasta(sequences = cseq, names = names(cseq), nbchar = 60, file.out = foc)
  pseq = proseq[ids]
  fop = sprintf("%s/01.pro.seq/%s.fas", dirw, fam)
  write.fasta(sequences = pseq, names = names(pseq), nbchar = 60, file.out = fop)

  system(sprintf("clustalo -i 01.pro.seq/%s.fas -o 03.pro.aln/%s.aln --outfmt=clu --force --full --full-iter", fam, fam))
  system(sprintf("clustalw2 -infile=03.pro.aln/%s.aln -tree -outputtree=phylip -clustering=NJ -kimura", fam))
}

### create gene pairs (remote / tandem)
tc = data.frame()
for (fam in fams) {
  f_tree = sprintf("%s/03.pro.aln/%s.ph", dirw, fam)
  tree = read.tree(f_tree)
  ids = tree$tip.label
  
  tcs = data.frame(aid = ids[-length(ids)], bid = ids[-1], alen = NA, blen = NA,
    mat = NA, mis = NA, idty = NA, gcov = NA, stringsAsFactors = F)
  for (i in 1:nrow(tcs)) {
    aid = tcs$aid[i]; bid = tcs$bid[i]
    aseq = proseq[[aid]]; bseq = proseq[[bid]]
    
    qlen = as.numeric(nchar(aseq))
    tlen = as.numeric(nchar(bseq))
    pw = pairwiseAlignment(AAString(aseq), AAString(bseq), type = "global",
    substitutionMatrix = "BLOSUM62", gapOpening = -3, gapExtension = -1)
    mat = nmatch(pw)
    mis = nmismatch(pw)
    indel = nindel(pw)
    qgap = indel@insertion[2]
    tgap =indel@deletion[2]
    alen = nchar(pw)
    stopifnot(alen == mis + mat + qgap + tgap)
    qres = qlen - (mat + mis + tgap)
    tres = tlen - (mat + mis + qgap)
    qgap = qgap + tres
    tgap = tgap + qres
    len = alen + qres + tres
    stopifnot(len == mis + mat + qgap + tgap)
    
    tcs$alen[i] = qlen; tcs$blen[i] = tlen; tcs$mat[i] = mat; tcs$mis[i] = mis;
    tcs$idty[i] = mat / (mat + mis)
    tcs$gcov[i] = (qgap + tgap) / len
  }
  tc = rbind(tc, tcs)
}

tc = cbind(tc, remote = NA, stringsAsFactors = F)
for (i in 1:nrow(tc)) {
  aid = tc$aid[i]; bid = tc$bid[i]
  aidx = tg$idx[tg$id == aid]
  bidx = tg$idx[tg$id == bid]
  tc$remote[i] = !(abs(aidx - bidx) <= 2)
}

tp = tc[tc$idty >= 0.6 & tc$gcov < 0.3,]
table(tp$remote)

### calculate Ka/Ks
pids = as.vector(rbind(tp$aid, tp$bid))

foc = sprintf("%s/11.pair.cds.fas", dirw, fam)
write.fasta(sequences = cdsseq[pids], names = pids, nbchar = 60, file.out = foc)
fop = sprintf("%s/11.pair.pro.fas", dirw, fam)
write.fasta(sequences = proseq[pids], names = pids, nbchar = 60, file.out = fop)

system("python $git/bio-pipeline/synonymous_calculation/synonymous_calc.py 11.pair.pro.fas 11.pair.cds.fas > 12.kaks")

tk = read.table("12.kaks", header = T, sep = ",", as.is = T)

#write.table(cbind(tp, tk[-1]), file = '13.tbl', col.names = T, row.names = F, sep = "\t", quote = F)

ta = data.frame(remote = tp$remote, dS = tk$dS.ng, dN = tk$dN.ng, stringsAsFactors = F)
#ta = data.frame(remote = tp$remote, dS = tk$dS.yn, dN = tk$dN.yn, stringsAsFactors = F)
ta = ta[ta$dS > 0 & ta$dS <= 1,]
ta = cbind(ta, dS.bin = cut(ta$dS, breaks=c(0,0.25,0.5,1)), stringsAsFactors = F)

to = ddply(ta, .(dS.bin, remote), summarise, n = length(dS), dN.q50 = median(dN/dS), dN.q25 = quantile(dN/dS, 0.25), dN.q75 = quantile(dN/dS, 0.75))
tt = ddply(to, .(dS.bin), summarise, y = min(dN.q25), label = paste(n, collapse = " / "))

p1 = ggplot(to) + 
  geom_crossbar(aes(x = dS.bin, y = dN.q50, ymin = dN.q25, ymax = dN.q75, fill = remote),
    stat = 'identity', position = 'dodge', geom_params = list(width = 0.5)) + 
  scale_fill_brewer(name = "", palette = "Paired", labels = c("Tandem (local)", "Ectopic")) +
  scale_x_discrete(name = 'dS bins') +
  scale_y_continuous(name = 'dN/dS') +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank(), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(1, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_text(size = 9, angle = 0), legend.text = element_text(size = 9, angle = 0)) +
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 90, hjust = 1)) +
  annotate("text", x = tt$dS.bin, y = tt$y, label = tt$label, vjust = 1.2, size = 3)

fo = sprintf("%s/12.%s.kaks.pdf", dirw, org)
ggsave(p1, filename = fo, width = 5, height = 5)



p1 = ggplot(ta, aes(x = dS, y = dN, shape = remote, color = remote)) + 
  geom_point() + 
  stat_smooth(method = lm, alpha = 0.2) +
  scale_color_manual(name = "Remote", values = c("forestgreen", "salmon"), labels = c("Tandem (local)", "Ectopic")) +
  scale_shape_manual(name = "Remote", values = c(19, 17), labels = c("Tandem (local)", "Ectopic")) +
  scale_x_continuous(name = 'dS') +
  scale_y_continuous(name = 'dN') +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(legend.position = c(0.8, 0.2), legend.title = element_blank(), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(1, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_text(size = 9, angle = 0), legend.text = element_text(size = 9, angle = 0)) +
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 90, hjust = 1))

fo = sprintf("%s/13.%s.kaks.pdf", dirw, org)
ggsave(p1, filename = fo, width = 5, height = 5)


## assign clades (obsolete)
#  cls = prosperi.cluster(tree, 0.2)$membership
#  ids = tree$tip.label

#  ff = sprintf("%s/04.cluster/%s.png", dirw, fam)
  
#  ucl = sort(unique(cls))
#  pal = rainbow(length(ucl))
#  names(pal) = as.character(ucl)

#  ht = length(tree$tip.label) * 12
#  png(filename = ff, width = 800, height = ht, units = 'px')
#  plot(tree, font = 1, show.node.label = F, show.tip.label = F,
#    label.offset = 0.01, no.margin = T, cex = 0.81)
#  newtip = paste(cls, tree$tip.label, sep = ' ')
#  tiplabels(newtip, adj = -0.1, bg = pal[as.character(cls)], cex = 0.8)

#  add.scale.bar(lcol = 'black')
#  dev.off()
  
  tcs = data.frame(fam = fam, id = ids, cl = cls, stringsAsFactors = F)
  tc = rbind(tc, tcs)
}

tc2 = tc[tc$cl > 0, ]
tc3 = merge(tc2, tg2[,c('id','chr','beg','end')], by = c('id'))

tp = data.frame()
for (fam in unique(tc3$fam)) {
  tcs = tc3[tc3$fam == fam,]
  stopifnot(nrow(tcs) >= 2)
  for (i in 2:nrow(tcs)) {
    id1 = tcs$id[i-1]; id2 = tcs$id[i]; cl = tcs$cl[i]
    remote = ifelse(
      tcs$chr[i-1] == tcs$chr[i] & abs(tcs$beg[i-1] - tcs$beg[i]) <= 200000,
      0, 1)
    tps = data.frame(fam = fam, id1 = id1, id2 = id2, cl = cl, rem = remote, stringsAsFactors = F)
    tp = rbind(tp, tps)
  }
}
tps = tp[tp$fam >= "CRP1020" & tp$fam < "CRP1600",]
table(tps$rem)

fo = sprintf("%s/27.cl.tbl", dirw)
write.table(to, file = fo, sep = "\t", row.names = F, col.names = T, quote = F)



clustalo -i 01.pro.seq/CRP0000.fas -o 03.pro.aln/CRP0000.aln --outfmt=clu --force --full --full-iter
pal2nal.pl 03.pro.aln/CRP0000.aln 01.cds.seq/CRP0000.fas -output paml > 03.cds.aln/CRP0000.phy
pal2nal.pl 03.pro.aln/CRP0000.aln 01.cds.seq/CRP0000.fas -output fasta > 03.cds.aln/CRP0000.fas

clustalw2 -INFILE=03.pro.aln/CRP0000.aln -BOOTSTRAP=1000 -OUTORDER=INPUT -OUTPUTTREE=newick -BOOTLABELS=node -CLUSTERING=NJ -KIMURA

    N-G Ka Ks
    Y-N Ka Ks