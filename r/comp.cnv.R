require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(rtracklayer)
require(hash)
require(grid)
source('Location.R')
source('comp.fun.R')

diro = file.path(Sys.getenv("misc3"), "comp.stat")

### looking at C.V. of family members
tg = data.frame()
for (qname in qnames_15) {
  cfg = cfgs[[qname]]
  fg = file.path(cfg$gdir, "51.gtb")
  tg1 = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1,16,17)]
  colnames(tg1)[2:3] = c('fam', 'subfam')
  tg = rbind(tg, cbind(org = qname, tg1))
}

to = tg
fams = c("CC-NBS-LRR", "TIR-NBS-LRR", "F-box", "LRR-RLK", "NCR", "TE", "Unknown")
fams = c(fams, "CRP0000-1030", "CRP1600-6250", "RLK")
to$fam[! to$fam %in% fams] = 'Pfam-Miscellaneous'
to$fam = factor(to$fam, levels = c(fams, 'Pfam-Miscellaneous'))

gb = group_by(to, org, fam, subfam)
to1 = summarise(gb, cnt = n())
to2 = ddply(to1, .(fam, subfam), summarise, cv = sd(cnt) / mean(cnt))
to3 = ddply(to2, .(fam), summarise, cnts = length(subfam), q25 = quantile(cv, 0.25, na.rm = T), q50 = median(cv, na.rm = T), q75 = quantile(cv, 0.75, na.rm = T))

fams = to3$fam[order(to3$q50, decreasing = T)]
to3$fam = factor(to3$fam, levels = fams)
labs = sprintf("%s (%d)", fams, to3$cnts[match(fams, to3$fam)])

p1 = ggplot(to3) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', geom_params = list(width = 0.5)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01), breaks = fams, labels = labs) +
  scale_y_continuous(name = 'C.V. of sub-family size') +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  theme(plot.margin = unit(c(0.5,1,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8))
  
fp = file.path(diro, "44.genefam.cv.pdf")
ggsave(p1, filename = fp, width = 5, height = 3)

### find syntenic ortholog
fg = file.path(Sys.getenv("genome"), tname, "51.tbl")
tg = read.table(fg, header = F, sep = "\t", as.is = T)
colnames(tg) = c('chr','beg','end','srd','id','type','fam')
tg = tg[tg$type == 'cds', c('chr','beg','end','id','fam')]
grt = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gb = group_by(tg, id)
tgl = summarise(gb, len = sum(end - beg + 1))

qname = "HM058"

fg = file.path(Sys.getenv("genome"), qname, "51.tbl")
qg = read.table(fg, header = F, sep = "\t", as.is = T)
colnames(qg) = c('chr','beg','end','srd','id','type','fam')
qg = qg[qg$type == 'cds', c('chr','beg','end','id','fam')]
grq = with(qg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gb = group_by(qg, id)
qgl = summarise(gb, len = sum(end - beg + 1))

fchain = sprintf("%s/%s_HM101/23_blat/31.9.chain", Sys.getenv("misc3"), qname)
chain = import.chain(fchain)

#x = tg[tg$id == 'Medtr1g004990.1',]
#grt = with(x, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
y = liftOver(grt, chain)
y = as(y, 'data.frame')

##### get cluster size for all accessions
diri = file.path(Sys.getenv("misc2"), "gene.cluster")

fc = sprintf("%s/11.tandem/%s.tbl", diri, tname)
tc = read.table(fc, header = T, sep = "\t", as.is = T)
tc = tc[,c('id','chr','beg','end','cat2','cat3','clu')]
colnames(tc)[5:6] = c('fam','sfam')
idxs = which(is.na(tc$clu))
tc$clu[idxs] = seq(max(tc$clu, na.rm = T)+1, by = 1, length.out = length(idxs))

gb = group_by(tc, clu)
tcl = summarise(gb, sfam = names(sort(table(sfam), decreasing = T))[1], fam = fam[which(sfam == sfam)[1]], chr = chr[1], beg = min(beg), end = max(end), csize = n(), nfam = length(unique(fam)))
grc = with(tcl, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

qname = 'HM058'
for (qname in qnames_all) {
cat(qname, "\n")
fa = sprintf("%s/%s_HM101/23_blat/31.9.gal", Sys.getenv("misc3"), qname)
ta = read.table(fa, header = T, sep = "\t", as.is = T)[,c('tId','tBeg','tEnd')]
gra = with(ta, GRanges(seqnames = tId, ranges = IRanges(tBeg, end = tEnd)))
gra = reduce(gra)

fc = sprintf("%s/11.tandem/%s.tbl", diri, qname)
qc = read.table(fc, header = T, sep = "\t", as.is = T)
qc = qc[,c('id','chr','beg','end','cat2','cat3','clu')]
colnames(qc)[5:6] = c('fam','sfam')
idxs = which(is.na(qc$clu))
qc$clu[idxs] = seq(max(qc$clu, na.rm = T)+1, by = 1, length.out = length(idxs))

gb = group_by(qc, clu)
qcl = summarise(gb, sfam = names(sort(table(sfam), decreasing = T))[1], fam = fam[which(sfam == sfam)[1]], chr = chr[1], beg = min(beg), end = max(end), csize = n(), nfam = length(unique(fam)))


#fchain = sprintf("%s/%s_HM101/23_blat/31.9.chain", Sys.getenv("misc3"), qname)
#chain = import.chain(fchain)

fh = sprintf("%s/%s_HM101/51_ortho/01.ortho.tbl", Sys.getenv("misc3"), qname)
th = read.table(fh, header = T, sep = "\t", as.is = T)
th = th[th$qid != '' & th$slen / pmin(th$qlen, th$tlen) >= 0.5,]
ths = th[,c('tid','qid')]
#h1 = hash(keys = th$tid, values = th$qid)
#h2 = hash(keys = th$qid, values = th$tid)

olens = intersect_basepair(grc, gra)
clus = tcl$clu[tcl$csize > 1 & olens/(tcl$end-tcl$beg+1) >= 0.95]
t1 = merge(tc[tc$clu %in% clus,], ths, by.x = 'id', by.y = 'tid')
t2 = merge(t1, qc, by.x = 'qid', by.y = 'id')
t3 = unique(t2[,c('clu.x','clu.y')])
colnames(t3) = c('tclu', 'qclu')

t4 = merge(t3, qc[,c('id','clu')], by.x = 'qclu', by.y = 'clu')
colnames(t4)[3] = 'qid'
t5 = merge(t4, ths, by = 'qid', all.x = T)
t6 = merge(t5, tc[,c('id','clu')], by.x = 'tid', by.y = 'id', all.x = T)
colnames(t6)[5] = 'tclu2'
t7 = cbind(t6, isSyn = is.na(t6$tid) | t6$tclu == t6$tclu2 | !duplicated(t6$qid))

gb = group_by(t7, tclu, qclu)
d1 = summarise(gb, nq = sum(isSyn))
d2 = merge(d1, tcl, by.x = 'tclu', by.y = 'clu')
gb = group_by(d2, tclu)
d3 = summarise(gb, nt = unique(csize), nq = sum(nq), fam = fam[1])

#d4 = cbind(d3, tag = d3$nt != d3$nq)
#ddply(d4, .(fam), summarise, cnt = sum(tag))

fo = sprintf("%s/31.cnv/%s.tbl", diri, qname)
write.table(d3[,c('tclu','nq')], fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')
}


##### merge cluster size and calculate CV
diri = file.path(Sys.getenv("misc2"), "gene.cluster")

fc = sprintf("%s/11.tandem/%s.tbl", diri, tname)
tc = read.table(fc, header = T, sep = "\t", as.is = T)
tc = tc[,c('id','chr','beg','end','cat2','cat3','clu')]
colnames(tc)[5:6] = c('fam','sfam')
idxs = which(is.na(tc$clu))
tc$clu[idxs] = seq(max(tc$clu, na.rm = T)+1, by = 1, length.out = length(idxs))

gb = group_by(tc, clu)
tl = summarise(gb, sfam = names(sort(table(sfam), decreasing = T))[1], fam = fam[which(sfam == sfam)[1]], HM101 = n())
tl = tl[tl$HM101 > 2,]

for (qname in qnames_12) {
  fc = sprintf("%s/31.cnv/%s.tbl", diri, qname)
  tc = read.table(fc, header = T, sep = "\t", as.is = T)
  tl = merge(tl, tc, by.x = 'clu', by.y = 'tclu', all.x = T)
  colnames(tl)[ncol(tl)] = qname
}

norg = apply(tl[,4:16], 1, cntorg <- function(x) sum(!is.na(x)))
to = tl[norg >= 8,]
cv = apply(to[,4:16], 1, cv <- function(x) sd(x, na.rm = T) / mean(x, na.rm = T))
to = cbind(to[,c('fam','sfam')], cv = cv)

fo = file.path(diro, "48.cnv.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

tp = ddply(to, .(fam), summarise, cnts = length(fam), q25 = quantile(cv, 0.25, na.rm = T), q50 = median(cv, na.rm = T), q75 = quantile(cv, 0.75, na.rm = T))

fams = tp$fam[order(tp$q50, decreasing = T)]
tp$fam = factor(tp$fam, levels = fams)
labs = sprintf("%s (%d)", fams, tp$cnts[match(fams, tp$fam)])

p1 = ggplot(tp) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', geom_params = list(width = 0.5)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01), labels = labs) +
  scale_y_continuous(name = 'C.V. of cluster size') +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  theme(plot.margin = unit(c(1,1,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8))
  
fp = file.path(diro, "48.cnv.pdf")
ggsave(p1, filename = fp, width = 5, height = 8)
