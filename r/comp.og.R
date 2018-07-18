require(plyr)
require(dplyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(reshape2)
require(RColorBrewer)
require(ape)
require(gtable)
require(gridBase)
require(colorRamps)
require(utils)
source("Location.R")
source("ggplot.fun.R")
source("comp.fun.R")

dirw = file.path('/home/youngn/zhoux379/data/misc3', "comp.og")
diro = file.path('/home/youngn/zhoux379/data/misc3', "comp.genefam")

## read all gene models
tg = data.frame()
for (qname in c(tname, qnames_12)) {
  gdir = cfgs[[qname]]$gdir
  fg = file.path(gdir, "51.gtb")
  ti = read.table(fg, sep = "\t", header = T, as.is = T)
  tg = rbind(tg, data.frame(org = qname, gid = ti$id, fam = ti$cat2, sfam = ti$cat3, stringsAsFactors = F))
}
table(tg$org)
fo = file.path(dirw, "01.gid.tbl")
#write.table(tg, fo, sep = "\t", row.names = F, col.names = T, quote = F)

## read gene cluster information
fi = file.path(dirw, "05.clu.70/32.tbl")
fi = file.path(dirw, "09.fastortho/11.tbl")
ti = read.table(fi, sep = "\t", header = T, as.is = T)
x = strsplit(ti$id, "[|]")
t_clu = cbind(ti, org = sapply(x, "[", 1), gid = sapply(x, "[", 2))

##### pan-proteome splited by gene-family
ti = t_clu
ti = ti[ti$org %in% c(tname,qnames_12),]
ti = merge(ti, tg, by = c('org', 'gid'))
ti = ti[order(ti$grp, ti$org), c(3,4,1,2,5,6)]
fo = file.path(dirw, "06.clu.tsv")
#write.table(ti, fo, sep = "\t", row.names = F, col.names = T, quote = F)

#ti = ti[ti$fam!='TE',]
gb = group_by(ti, grp)
tr = dplyr::summarise(gb, size = length(unique(org)), org = org[1], fam = names(sort(table(fam), decreasing = T))[1], rid = id[1])
ti = merge(ti, tr[,c('grp','size')], by = 'grp')

# output unknown cluster/ids for blastnr
ids_unk = tr$rid[tr$fam == 'Unknown']
fo = file.path(dirw, "31.unk.txt")
#write(ids_unk, fo, sep = "\n")
#seqret.pl -d 02.fas -b 31.unk.txt -o 32.unk.fas

table(tr$size)
table(tr$org[tr$size==1])

grps = tr$grp[tr$size == 1]
ts = tr[tr$grp %in% grps,]

y = table(ts$fam)
y[order(y, decreasing = T)][1:20]/sum(y)

brks = as.character(1:13)
# plot all fams -hist

tq = ddply(tr, .(fam), summarise, cnt = length(fam), avg = mean(size))
fams = tq$fam[order(tq$avg, decreasing=T)]

do = data.frame()
for (fam in fams) {
  x = table(tr$size[tr$fam == fam])
  dos = data.frame(fam = fam, size = brks, cnt = as.numeric(x[brks]), stringsAsFactors = F)
  do = rbind(do, dos)
}
do$cnt[is.na(do$cnt)] = 0
do$size = factor(do$size, levels = brks)
do$fam = factor(do$fam, levels = fams)

p1 = ggplot(do, aes(x = size, y = cnt)) +
  geom_bar(stat = 'identity', width = 0.8) +
  scale_x_discrete(name = '', breaks = brks) +
  scale_y_continuous(name = '# Clusters') +
  facet_wrap(~ fam, scales = 'free', nrow = 10) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 90, hjust = 1))

fp = sprintf("%s/22.afs.pdf", diro)
ggsave(p1, filename = fp, width = 16, height = 20)

# plot select fams -hist
to = rename_genefam(tr)
tq = ddply(to, .(fam), summarise, cnt = length(fam), avg = mean(size))
fams = tq$fam[order(tq$avg, decreasing=T)]

do = data.frame()
for (fam in fams) {
  x = table(to$size[to$fam == fam])
  dos = data.frame(fam = fam, size = brks, cnt = as.numeric(x[brks]), stringsAsFactors = F)
  do = rbind(do, dos)
}
do$cnt[is.na(do$cnt)] = 0
do$size = factor(do$size, levels = brks)
do$fam = factor(do$fam, levels = fams)

p1 = ggplot(do, aes(x = size, y = cnt)) +
  geom_bar(stat = 'identity', width = 0.8) +
  scale_x_discrete(name = '', breaks = brks) +
  scale_y_continuous(name = '# Clusters') +
  facet_wrap(~ fam, scales = 'free', nrow = 4) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 90, hjust = 1))

fp = sprintf("%s/22.afs.select.pdf", diro)
ggsave(p1, filename = fp, width = 6, height = 8)

dos = ddply(do, .(fam), summarise, total = sum(cnt))
to = do[do$size %in% c(1,13),]
to$size = factor(to$size, levels = unique(to$size))
to = merge(to, dos, by = 'fam')
to = cbind(to, prop = to$cnt/to$total)
tp = ddply(to, .(fam), summarise, fold = prop[which(size==1)]/prop[which(size==13)])
fams = tp$fam[order(tp$fold)]
tp$fam = factor(tp$fam, levels = fams)
to$fam = factor(to$fam, levels = fams)

fo = file.path(dirw, "41.bar.select.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

cols = brewer.pal(3,"Dark2")[1:2]
p1 = ggplot(to, aes(x = fam, y = prop, fill = size)) +
  geom_bar(stat = 'identity', position = 'dodge', width = 0.6) +
  scale_x_discrete(name = '') +
  scale_y_continuous(name = 'Proportion of gene family', expand=c(0,0), limits=c(0,0.8)) +
  scale_fill_manual(labels=c("Shared by only 1 accession (accession-specific)", "Shared by all 13 accessions"), values=cols) +
  theme_bw() +
  theme(legend.position=c(0.5,0.9), legend.background = element_rect(fill=NA, colour=NA, size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.6, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_blank(), legend.text = element_text(size = 9, angle = 0)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.8,0.5,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size=8, color="black", angle=0, hjust=0.5))

p2 = ggplot(tp, aes(x = fam, y = fold)) +
  geom_point(size = 0.9) +
  geom_line(aes(group=1), size = 0.3) +
  scale_x_discrete(name = '') +
  scale_y_continuous(name = 'Fold change', expand=c(0,0), limits=c(0,6.2)) +
  theme_bw() + 
#  theme(panel.grid = element_blank()) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0,0.5,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size=8, color="blue", angle=30, hjust=1)) +
  theme(axis.text.y = element_text(size=8, color="black", angle=0, hjust=0.5))

g1 = ggplotGrob(p1)
g2 = ggplotGrob(p2)
g2$widths = g1$widths
gs = list(g1, g2)
heis = c(4,3)
g <- gtable_matrix(name = 'demo', grobs = matrix(gs, nrow = length(gs)), widths = 1, heights = heis)

fo = sprintf("%s/22.bar.select.pdf", diro)
pdf(file = fo, width = 5, height = 4, bg = 'transparent')
grid.newpage()
grid.draw(g)
dev.off()


# select fams - violin
to = rename_genefam(tr)
tos = ddply(to, .(fam), summarise, size_median = median(size))
fams = tos$fam[order(tos$size_median, decreasing = T)]
to$fam = factor(to$fam, levels = fams)

p1 = ggplot(to) +
  geom_violin(aes(x = fam, y = size), draw_quantiles = c(0.25, 0.5, 0.75)) + 
  coord_flip() +
  scale_x_discrete(name = '') +
  scale_y_continuous(name = '# Shared Accession', expand = c(0, 0), limits = c(0.5, 13.5)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))
fp = sprintf("%s/22.afs.select.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 8)


# CNV index
gb = group_by(ti[ti$size>1,], grp, org)
tv = dplyr::summarise(gb, cn = length(org))
t1 = expand.grid(grp=unique(ti$grp), org=unique(ti$org))
t2 = merge(t1, tv, by=c('grp','org'), all.x=T)
t2$cn[is.na(t2$cn)]=0
gb = group_by(t2, grp)
to = dplyr::summarise(gb, cv = sd(cn) / mean(cn))
to = merge(to, tr[,c('grp','fam')], by = 'grp')

fo = file.path(diro, "25.cnv.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

tp = ddply(to, .(fam), summarise, cnt = length(fam), q25=quantile(cv,0.25,na.rm=T), q50=quantile(cv,0.5,na.rm=T), q75=quantile(cv,0.75,na.rm=T), avg=mean(cv,na.rm=T), std=sd(cv,na.rm=T))
tp = tp[order(tp$q50, tp$q75, decreasing=T),]
fams = tp$fam
tp$fam = factor(tp$fam, levels = fams)

p1 = ggplot(tp) +
	#geom_errorbar(aes(x=fam,ymin=q25,ymax=q75), width=0.2) +
  #geom_point(aes(x=fam,y=q50), size=0.7) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', position = 'dodge', width = 0.6) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01), labels = sprintf("%s | %5d", tp$fam, tp$cnt)) +
  scale_y_continuous(name = 'CNV Index', expand = c(0.01, 0.01)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0,'lines')) +
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = file.path(diro, "25.cnv.pdf")
ggsave(p1, filename = fp, width = 5, height = 8)


##### enrichment of unknown proteins in low-confidence calls
fq = file.path(Sys.getenv("genome"), "HM101", "augustus", "31.gtb")
tq = read.table(fq, sep = "\t", header = T, as.is = T)

idxs = which(grepl("\\[\\w+\\]", tq$note))
m = regexpr("\\[(\\w+)\\]", tq$note, perl = TRUE)
x = regmatches(tq$note, m)
stopifnot(length(x) == length(idxs))
tqs = data.frame(id = tq$id[idxs], qual = x, note = tq$note[idxs])
ids_hc = paste(tname, tqs$id[tqs$qual == '[HC]'], sep = "|")
ids_lc = paste(tname, tqs$id[tqs$qual == '[LC]'], sep = "|")

x = strsplit(tr$rid, "|")
gids = sapply(x, "[", 2)
#table(tqs$qual[tqs$id %in% gids[tr$org == 'HM101' & tr$size >1 & tr$fam != 'Unknown']])


ds = data.frame()
for (qname in qnames_12) {
  fq = file.path(Sys.getenv("genome"), qname, "augustus", "31.gtb")
  tq = read.table(fq, sep = "\t", header = T, as.is = T)
  dss = tq[,c('id','note')]
  colnames(dss)[2] = 'qual'
  dss = cbind(dss, rid = paste(qname, dss$id, sep = "|"), stringsAsFactors = F)
  ds = rbind(ds, dss)
}
#hist(ds$qual[ds$rid %in% ti$id[ti$size == 1 & ti$cat2 == 'Unknown']])
#hist(ds$qual[ds$rid %in% ti$id[ti$size == 1 & ti$cat2 != 'Unknown']])
#hist(ds$qual[ds$rid %in% ti$id[ti$size == 13 & ti$cat2 == 'Unknown']])
ids_hc = c(ids_hc, ds$rid[ds$qual >= 0.9])
ids_lc = c(ids_lc, ds$rid[ds$qual < 0.9])


do = data.frame()
for (qname in c("HM056")) {
  fx = file.path(Sys.getenv("misc2"), "rnaseq/mt/31_cufflinks", qname, "isoforms.fpkm_tracking")
  tx = read.table(fx, sep = "\t", header = T, as.is = T)
  dos = data.frame(rid = paste(qname, tx$tracking_id, sep = "-"), fpkm = tx$FPKM, stringsAsFactors = F)
  do = rbind(do, dos)
}
do$fpkm[do$fpkm > 0] = 1
#hist(do$fpkm[do$rid %in% ti$id[ti$size == 1 & ti$cat2 != 'Unknown']])
#hist(do$fpkm[do$rid %in% ti$id[ti$size == 1 & ti$cat2 == 'Unknown']])
#hist(do$fpkm[do$rid %in% ti$id[ti$size > 12 & ti$cat2 == 'Unknown']])


t.unk = tr[tr$fam == 'Pfam:other',]
t.unk = tr[tr$fam == 'Unknown',]
t.unk.hc = t.unk[t.unk$rid %in% ids_hc,]
t.unk.lc = t.unk[t.unk$rid %in% ids_lc,]

brks = as.character(1:13)

x = table(t.unk$size)
do1 = data.frame(lab = 'All Unknown Proteins', size = brks, cnt = as.numeric(x[brks]), stringsAsFactors = F)
x = table(t.unk.hc$size)
do2 = data.frame(lab = 'High Confidence', size = brks, cnt = as.numeric(x[brks]), stringsAsFactors = F)
x = table(t.unk.lc$size)
do3 = data.frame(lab = 'Low Confidence', size = brks, cnt = as.numeric(x[brks]), stringsAsFactors = F)
do = rbind(do2, do3)
do$size = factor(do$size, levels = brks)

p1 = ggplot(do, aes(x = size, y = cnt, fill = lab)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.8) +
  scale_x_discrete(name = '', breaks = brks) +
  scale_y_continuous(name = '# Clusters') +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = c(0.6, 0.8), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.6, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_blank(), legend.text = element_text(size = 9, angle = 0)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 90, hjust = 1))

fp = sprintf("%s/45.unk.pdf", dirw)
ggsave(p1, filename = fp, width = 4, height = 4)


##### accession-specific / core genes - pan-proteome approach
ti = t_clu
ti = ti[ti$org %in% c(tname,qnames_12),]
ti = merge(ti, tg, by = c('org', 'gid'))
#ids = sprintf("%s-%s", ti$org, ti$gid)
#ti = ti[! ids %in% ids_low,]

gb = group_by(ti, grp)
tr = dplyr::summarise(gb, size = length(unique(org)), org = org[1], fam = names(sort(table(fam), decreasing = T))[1], rid = id[1])
trs = tr[tr$fam != 'TE',]

tas = trs[trs$size==1,]
tco = trs[trs$size==13,]
table(tas$fam)/nrow(tas)
table(tco$fam)/nrow(tco)

fams = as.character(unique(trs$fam))

x1 = table(tas$fam)
num1 = as.numeric(x1[fams])
tp1 = data.frame(opt = 'pan-proteome', type = 'accession-specific', fam = fams, num = num1, total = nrow(tas), stringsAsFactors = T)

x2 = table(tco$fam)
num2 = as.numeric(x2[fams])
tp2 = data.frame(opt = 'pan-proteome', type = 'core', fam = fams, num = num2, total = nrow(tco), stringsAsFactors = T)

##### accession-specific / core genes - pan-genome approach
## core-genome genes + HM101-specific genes
tlen = read.table(cfgs[[tname]]$size, sep = "\t", as.is = T, header = F)
tgap = read.table(cfgs[[tname]]$gap, sep = "\t", as.is = T, header = F)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, V2)))
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, V3)))
grs = setdiff(grt, grp)

gr_core = grs; gr_as = grs
for (org in qnames_12) {
  fgax = sprintf("%s/%s_HM101/23_blat/31.9/gax", Sys.getenv('misc3'), org)
  tgax = read.table(fgax, sep = "\t", header = F, as.is = T)[,1:3]
  gr = with(tgax, GRanges(seqnames = V1, ranges = IRanges(V2, V3)))
  gr_core = intersect(gr_core, gr)
  gr_as = setdiff(gr_as, gr)
}

  td = read.table(tcfg$gene, sep = "\t", header = F, as.is = T)
  colnames(td) = c("chr", "beg", "end", "srd", "id", "type", "fam")
  td = td[td$type == 'cds',]
  grg = with(td, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

  olens = intersect_basepair(grg, gr_as)
  gb = group_by(cbind(td, olen = olens), id)
  dx = dplyr::summarise(gb, len = sum(end-beg+1), olen = sum(olen), pct = olen/len)
  dy = dx[dx$pct >= 0.5,]
t_asr = data.frame(org = tname, gid = dy$id)

  olens = intersect_basepair(grg, gr_core)
  gb = group_by(cbind(td, olen = olens), id)
  dx = dplyr::summarise(gb, len = sum(end-beg+1), olen = sum(olen), pct = olen/len)
  dy = dx[dx$pct >= 0.5,]
t_core = data.frame(org = tname, gid = dy$id)

## acc12-specific genes
fi = file.path(Sys.getenv('misc3'), 'comp.panseq', '32.global.tbl')
t_aln = read.table(fi, header = T, sep = "\t", as.is = T)

gb = group_by(t_aln, cid)
dcl = dplyr::summarise(gb, n_org = n(), size = sum(end - beg + 1))
cids = dcl$cid[dcl$n_org == 1]

do = data.frame()
for (qname in qnames_12) {
  cfg = cfgs[[qname]]

  td = read.table(cfg$gene, sep = "\t", header = F, as.is = T)
  colnames(td) = c("chr", "beg", "end", "srd", "id", "type", "fam")
  td = td[td$type == 'cds',]
  grg = with(td, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

  fn = file.path(cfg$cdir, "../41_novseq/21.bed")
#  tn = read.table(fn, sep = "\t", header = F, as.is = T)
#  grn = with(tn, GRanges(seqnames = V1, ranges = IRanges(V2+1, end = V3)))
  tn = t_aln[t_aln$cid %in% cids & t_aln$org == qname,]
  grn = with(tn, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

  olens = intersect_basepair(grg, grn)
  #ds = ddply(cbind(tg, olen = olens), .(fam), summarise, olen = sum(olen), alen = sum(end - beg + 1))
  #ds = cbind(ds, prop = ds$olen / ds$alen, org = qname)
  gb = group_by(cbind(td, olen = olens), id)
  dx = dplyr::summarise(gb, len = sum(end-beg+1), olen = sum(olen), pct = olen/len)

  dy = dx[dx$pct >= 0.5,]
  do = rbind(do, cbind(dy, org = qname))
}
t_nov = do[,c('org','id')]
colnames(t_nov)[2] = 'gid'

tas = merge(t_nov, tg, by = c('org', 'gid'))
tas = tas[tas$fam != 'TE',]
table(tas$fam)/nrow(tas)

t2 = merge(t_clu, t_core, by = c('org','gid'))
t3 = merge(t2, tg, by = c('org', 'gid'))
gb = group_by(t3, grp)
tr = dplyr::summarise(gb, size = length(unique(org)), org = org[1], fam = names(sort(table(fam), decreasing = T))[1], rid = id[1])
tco = merge(t3, tr[,c('grp','size')], by = 'grp')
tco = tco[tco$fam != 'TE',]
table(tco$fam)/nrow(tco)


x1 = table(tas$fam)
num1 = as.numeric(x1[fams])
tg1 = data.frame(opt = 'pan-genome', type = 'accession-specific', fam = fams, num = num1, total = nrow(tas), stringsAsFactors = T)

x2 = table(tco$fam)
num2 = as.numeric(x2[fams])
tg2 = data.frame(opt = 'pan-genome', type = 'core', fam = fams, num = num2, total = nrow(tco), stringsAsFactors = T)

## save to file
to = rbind(tp1, tp2, tg1, tg2)
to = to[!is.na(to$num),]
fo = file.path(dirw, "61.core+nov.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)


### combo-piechart
fi = file.path(dirw, "61.core+nov.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)

split_pie_data <- function(tt) {
	bfams = c("Pfam:other", "Unknown")
	
	tt2 = tt[!tt$fam %in% bfams,]
	tt2$fam = factor(tt2$fam, levels = unique(as.character(tt2$fam)))
	tt2 = tt2[order(tt2$fam),]
	tt2 = cbind(tt2, pct = tt2$num/tt2$total * 100)
	tmp = cumsum(tt2$pct)
	tmp = c(0, tmp[-length(tmp)])
	tt2 = cbind(tt2, lab = sprintf("%.01f", tt2$pct), x = tt2$pct/2 + tmp)
	
	tt1 = tt
	tt1$fam[!tt1$fam %in% bfams] = 'ESGF'
	tt1$fam = factor(tt1$fam, levels = unique(as.character(tt1$fam)))
	tt1 = tt1[order(tt1$fam),]
	tt1 = ddply(tt1, .(fam), summarise, num=sum(num), total=total[1])
	tt1 = cbind(tt1, pct = tt1$num/tt1$total * 100)
	tmp = cumsum(tt1$pct)
	tmp = c(0, tmp[-length(tmp)])
	tt1 = cbind(tt1, lab = sprintf("%.01f", tt1$pct), x = tt1$pct/2 + tmp)
	
	list(tt1, tt2)
}
	
lst = split_pie_data(ti[ti$opt=='pan-genome' & ti$type=='core',])
tg1a = lst[[1]]; tg1b = lst[[2]]
lst = split_pie_data(ti[ti$opt=='pan-genome' & ti$type=='accession-specific',])
tg2a = lst[[1]]; tg2b = lst[[2]]
lst = split_pie_data(ti[ti$opt=='pan-proteome' & ti$type=='core',])
tp1a = lst[[1]]; tp1b = lst[[2]]
lst = split_pie_data(ti[ti$opt=='pan-proteome' & ti$type=='accession-specific',])
tp2a = lst[[1]]; tp2b = lst[[2]]


pg1a = ggplot(tg1a, aes(x = '', y = pct, fill = fam)) +
    geom_bar(stat = 'identity', width = 1) +
  	scale_y_continuous(name = '', breaks = tg1a$x, labels = tg1a$lab)
pg1b = ggplot(tg1b, aes(x = '', y = pct, fill = fam)) +
    geom_bar(stat = 'identity', width = 1) +
  	scale_y_continuous(name = '', breaks = tg1b$x, labels = tg1b$lab)
pg2a = ggplot(tg2a, aes(x = '', y = pct, fill = fam)) +
    geom_bar(stat = 'identity', width = 1) +
  	scale_y_continuous(name = '', breaks = tg2a$x, labels = tg2a$lab)
pg2b = ggplot(tg2b, aes(x = '', y = pct, fill = fam)) +
    geom_bar(stat = 'identity', width = 1) +
  	scale_y_continuous(name = '', breaks = tg2b$x, labels = tg2b$lab)
pp1a = ggplot(tp1a, aes(x = '', y = pct, fill = fam)) +
    geom_bar(stat = 'identity', width = 1) +
  	scale_y_continuous(name = '', breaks = tp1a$x, labels = tp1a$lab)
pp1b = ggplot(tp1b, aes(x = '', y = pct, fill = fam)) +
    geom_bar(stat = 'identity', width = 1) +
  	scale_y_continuous(name = '', breaks = tp1b$x, labels = tp1b$lab)
pp2a = ggplot(tp2a, aes(x = '', y = pct, fill = fam)) +
    geom_bar(stat = 'identity', width = 1) +
  	scale_y_continuous(name = '', breaks = tp2a$x, labels = tp2a$lab)
pp2b = ggplot(tp2b, aes(x = '', y = pct, fill = fam)) +
    geom_bar(stat = 'identity', width = 1) +
  	scale_y_continuous(name = '', breaks = tp2b$x, labels = tp2b$lab)


p_tmp = pg1b + 
  	scale_fill_manual(breaks = fams, values = cols) +
  	theme(legend.position = 'right', legend.background = element_rect(fill='NA',linetype=0), legend.key = element_rect(fill=NA, colour=NA, size=0), legend.key.size = unit(0.6,'lines'), legend.margin = unit(0,"lines"), legend.title = element_blank(), legend.text = element_text(size=9, angle=0))
leg = get_ggplot_legend(p_tmp)

famsa = levels(tg1a$fam)
famsb = levels(tg1b$fam)
colsa = c('grey', brewer.pal(3,"Set3")[1:2])
colsb = brewer.pal(12,"Set3")[3:12]
fams = c(famsa[c(2,3,1)],famsb)
cols = c(colsa[c(2,3,1)],colsb)

numrow = 4; hts = c(2,2,2,2)
numcol = 3; wds = c(2,1.5,1)
fp = sprintf("%s/61.core+nov.pdf", dirw)
pdf(file = fp, width = sum(wds), height = sum(hts), bg = 'transparent')
grid.newpage()
pushViewport(viewport(layout = grid.layout(numrow, numcol, width = wds, heights = hts)))

pas = list(pg1a,pg2a,pp1a,pp2a)
pbs = list(pg1b,pg2b,pp1b,pp2b)
xs = c(1,2,3,4)
ss = LETTERS[1:length(pas)]

for (i in 1:length(pas)) {
  x = xs[i]; lab = ss[i]
  pa = pas[[i]] +
  	scale_fill_manual(breaks = famsa, values = colsa) +
  	coord_polar(theta = "y") +
  	theme_bw() +
  	theme(panel.border = element_rect(fill=NA, linetype=0)) +
  	theme(legend.position = 'none') +
  	theme(axis.ticks.length = unit(0, 'lines')) +
  	theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  	theme(axis.title.x = element_blank()) +
  	theme(axis.title.y = element_blank()) +
  	theme(axis.text.x = element_text(size = 8, colour = "blue"))
  pb = pbs[[i]] +
  	scale_fill_manual(breaks = famsb, values = colsb) +
  	coord_polar(theta = "y") +
  	theme_bw() +
  	theme(panel.border = element_rect(fill=NA, linetype=0)) +
  	theme(legend.position = 'none') +
  	theme(axis.ticks.length = unit(0, 'lines')) +
  	theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  	theme(axis.title.x = element_blank()) +
  	theme(axis.title.y = element_blank()) +
  	theme(axis.text.x = element_text(size = 8, colour = "blue"))
  
  print(pa, vp = viewport(layout.pos.row = x, layout.pos.col = 1))
  print(pb, vp = viewport(layout.pos.row = x, layout.pos.col = 2))
  
  grid.text(lab, x = unit(0.1,'lines'), y = unit(1,'npc')-unit(0.1,'lines'), just = c('left', 'top'), gp = gpar(col = "black", fontface = 2, fontsize = 16),
    vp = viewport(layout.pos.row = x, layout.pos.col = 1))
}
grid.draw(legendGrob(famsa[c(2,3,1)], pch=15, hgap=unit(0.2,'lines'), vgap=unit(0.1,'lines'), gp=gpar(col=colsa[c(2,3,1)], fontsize=9), vp=viewport(layout.pos.row = c(2), layout.pos.col = 3)))
grid.draw(legendGrob(famsb, pch=15, hgap=unit(0.2,'lines'), vgap=unit(0.1,'lines'), gp=gpar(col=colsb, fontsize=9), vp=viewport(layout.pos.row = c(3), layout.pos.col = 3)))
dev.off()

## simple piechart
fi = file.path(dirw, "61.core+nov.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)

rename_genefam_nov <- function(ti) {
  mapping = list(
  	"Unknown" = "Unknown",
  	"ESGF" = c("CC-NBS-LRR", "TIR-NBS-LRR", "NB-ARC", "TIR", "CRP:NCR", "LRR", "HSP70", "Cytochrome", "Histone", 'Pkinase', 'Pkinase_Tyr', "RLK", "LRR-RLK")
  )
  to = ti
  for (fam in names(mapping)) {
  	to$fam[to$fam %in% mapping[[fam]]] = fam
  }
  to$fam[! to$fam %in% names(mapping)] = 'Pfam:other'
  to$fam = factor(to$fam, levels = unique(to$fam))
  to
}
ti = rename_genefam_nov(ti)
to = ddply(ti, .(opt, type, fam), summarise, num=sum(num), total=total[1])

to$fam = factor(to$fam, levels = unique(as.character(to$fam)))
to = to[order(to$fam),]
to = cbind(to, pct = to$num/to$total * 100)

do = ddply(to, .(opt, type), xt <- function(ds) {
  tmp = cumsum(ds$pct)
  tmp = c(0, tmp[-length(tmp)])
  cbind(ds, lab = sprintf("%.01f", ds$pct), x = ds$pct/2 + tmp)
})
#fo = file.path(dirw, "62.tbl")
#write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

p1 = ggplot(do, aes(x = '', y = pct, fill = fam)) +
  geom_bar(stat = 'identity', width = 1) +
  geom_text(data=do, mapping=aes(x=1.1, y=x, label=lab))+
  coord_polar(theta = "y") +
  facet_wrap(opt ~ type) +
#  scale_x_discrete(name = '', breaks = brks) +
  scale_y_continuous(name='') +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(legend.position = 'top', legend.background = element_rect(fill = 'white', colour=NA, size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.6, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_blank(), legend.text = element_text(size = 9, angle = 0)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank())

fo = sprintf("%s/43.nov.pdf", diro)
ggsave(p1, filename = fo, width = 5, height = 6)

