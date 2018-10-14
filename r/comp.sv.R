require(rtracklayer)
require(plyr)
require(grid)
require(ggplot2)
source('Location.R')
source('comp.fun.R')
source('comp.plot.fun.R')

diro = sprintf("%s/comp.stat", Sys.getenv("misc3"))

tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

grnp = GenomicRanges::setdiff(grt, grp)
tsize = sum(width(grt))
tsize2 = sum(width(grnp))

tg = read.table(tcfg$gene, sep = "\t", header = F, as.is = T)
colnames(tg) = c('chr', 'beg', 'end', 'srd', 'id', 'type', 'fam')
tg = tg[tg$type == 'cds',]
gb = group_by(tg, id)
tg = summarise(gb, fam = fam[1], chr = chr[1], beg = min(beg), end = max(end), size = end - beg + 1)
tg = cbind(idx = 1:nrow(tg), tg)
grgt = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

##### SV impact on genome and proteome
qname = "HM004"
for (qname in qnames_all) {
cfg = cfgs[[qname]]

tlen = read.table(cfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))

tgap = read.table(cfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

grnp = GenomicRanges::setdiff(grt, grp)
qsize = sum(width(grt))
qsize2 = sum(width(grnp))

qg = read.table(cfg$gene, sep = "\t", header = F, as.is = T)
colnames(qg) = c('chr', 'beg', 'end', 'srd', 'id', 'type', 'fam')
qg = qg[qg$type == 'cds',]
gb = group_by(qg, id)
qg = summarise(gb, fam = fam[1], chr = chr[1], beg = min(beg), end = max(end), size = end - beg + 1)
grgq = with(qg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

cdir = cfg$cdir
fy = file.path(cdir, "31.9/gal")
ty = read.table(fy, header = T, sep = "\t", as.is = T)[,c('tId', 'tBeg', 'tEnd', 'qId', 'qBeg', 'qEnd')]
gryt = with(ty, GRanges(seqnames = tId, ranges = IRanges(tBeg, end = tEnd)))
gryq = with(ty, GRanges(seqnames = qId, ranges = IRanges(qBeg, end = qEnd)))
ytlen = sum(width(reduce(gryt)))
yqlen = sum(width(reduce(gryq)))

fv = file.path(cdir, "../31_sv/05.stb")
tv = read.table(fv, header = T, sep = "\t", as.is = T)[,c('id','tchr','tbeg','tend','tlen','srd','qchr','qbeg','qend','qlen')]
tv = tv[tv$tlen+tv$qlen-2 >= 50,]
tvt = tv[tv$tlen > 0, c('tchr','tbeg','tend')]
grst = with(tvt, GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend)))
tvq = tv[tv$qlen > 0, c('qchr','qbeg','qend')]
grsq = with(tvq, GRanges(seqnames = qchr, ranges = IRanges(qbeg, end = qend)))

stlen = sum(width(reduce(grst)))
sqlen = sum(width(reduce(grsq)))

tocnt = intersect_count(grgt, grst)
tg2 = cbind(tg, cnt = tocnt)
tgb = group_by(tg2, id)
x = summarise(tgb, cnt = sum(cnt))
tpct = sum(x$cnt > 0) / nrow(x)

qocnt = intersect_count(grgq, grsq)
qg2 = cbind(qg, cnt = qocnt)
qgb = group_by(qg2, id)
x = summarise(qgb, cnt = sum(cnt))
qpct = sum(x$cnt > 0) / nrow(x)

cat(sprintf("%s: genome [tgt %.03f qry %.03f] gene [tgt %.03f qry %.03f]\n", qname, stlen/ytlen, sqlen/yqlen, tpct, qpct))

}

#### SV impact on genes / gene family
do = data.frame()
for (qname in qnames_all[1:12]) {

qname = "HM340.AC"
cfg = cfgs[[qname]]

tlen = read.table(cfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))

tgap = read.table(cfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

grnp = GenomicRanges::setdiff(grt, grp)
qsize = sum(width(grt))
qsize2 = sum(width(grnp))

qg = read.table(cfg$gene, sep = "\t", header = F, as.is = T)
colnames(qg) = c('chr', 'beg', 'end', 'srd', 'id', 'type', 'fam')
qg = qg[qg$type == 'cds',]
gb = group_by(qg, id)
qg = summarise(gb, fam = fam[1], chr = chr[1], beg = min(beg), end = max(end), size = end - beg + 1)
qg = cbind(idx = 1:nrow(qg), qg)
grgq = with(qg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

cdir = cfg$cdir
#fy = file.path(cdir, "31.9/idm")
#ty = read.table(fy, header = F, sep = "\t", as.is = T)
#colnames(ty) = c("tchr", 'tbeg', 'tend', 'tsrd', 'qchr', 'qbeg', 'qend', 'qsrd', 'cid', 'lev')
#ty = ty[ty$lev == 1,]
fv = file.path(cdir, "../31_sv/05.stb")
tv = read.table(fv, header = T, sep = "\t", as.is = T)[,c('id','tchr','tbeg','tend','tlen','srd','qchr','qbeg','qend','qlen')]
tv = tv[tv$tlen+tv$qlen-2 >= 50,]
tvt = tv[tv$tlen > 0, c('tchr','tbeg','tend')]
tvt = cbind(idx = 1:nrow(tvt), tvt)
grst = with(tvt, GRanges(seqnames = tchr, ranges = IRanges(tbeg+1, end = tend-1)))
tvq = tv[tv$qlen > 0, c('qchr','qbeg','qend','tchr','tbeg','tend')]
tvq = cbind(idx = 1:nrow(tvq), tvq)
grsq = with(tvq, GRanges(seqnames = qchr, ranges = IRanges(qbeg+1, end = qend-1)))

x = intersect_idx(grst, grgt)
y = merge(x, tg[,c('idx','id','fam','size')], by.x = 'qidx', by.y = 'idx')
gb = group_by(y, idx)
z = summarize(gb, ngene = n(), n2 = sum(fam %in% c("NCR","CRP0000-1030")))
z = merge(z, tvt, by = 'idx')
z[order(z$n2, decreasing = T), ][1:30,]

x = intersect_idx(grsq, grgq)
y = merge(x, qg[,c('idx','id','fam','size')], by.x = 'qidx', by.y = 'idx')
gb = group_by(y, idx)
z = summarize(gb, ngene = n(), n2 = sum(fam %in% c("NCR","CRP0000-1030")))
z = merge(z, tvq, by = 'idx')
z[order(z$n2, decreasing = T), ][1:30,]


ds = cbind(ds, org = qname)
do = rbind(do, ds)
}

to = ddply(do, .(fam), summarise, q25 = quantile(prop, 0.25), q50 = quantile(prop, 0.5), q75 = quantile(prop, 0.75))

ffam = file.path(diro, "41.gene.fams.tbl")
fams = read.table(ffam, header = F, sep = "\t", as.is = T)[,1]

ti = to
tis = ti[ti$fam %in% fams,]
tis = tis[order(tis$q50, decreasing = T),]
tis$fam = factor(tis$fam, levels = tis$fam)
p4 = ggplot(tis) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', position = 'dodge', geom_params = list(width = 0.7, size = 0.3)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01)) +
  scale_y_continuous(name = 'Proportion SV', expand = c(0, 0)) +
  theme_bw() +
#  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = file.path(dirw, "19_sv_genefam.pdf")
ggsave(p4, filename = fp, width = 5, height = 4)


