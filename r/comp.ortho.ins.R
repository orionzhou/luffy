require(plyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(RColorBrewer)
require(reshape2)
source("Location.R")
source("comp.fun.R")

dirw = file.path(Sys.getenv("misc3"), "comp.ortho.ins")

qname = 'HM004'
for (qname in qnames) {
qcfg = cfgs[[qname]]
ccfg = ccfgs[[qname]]

##### find ins/del genes
cfg = tcfg

tt = read.table(cfg$size, sep = "\t", header = F, as.is = T)
tt = data.frame(chr = tt$V1, beg = 1, end = tt$V2)
grt = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

tp = read.table(cfg$gap, sep = "\t", header = F, as.is = T)
colnames(tp) = c("chr", "beg", "end")
tp = tp[tp$end - tp$beg + 1 > 100,]
grp = with(tp, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

gru = GenomicRanges::setdiff(grt, grp)

###
fx = file.path(ccfg$cdir, "31.9/gax")
tx = read.table(fx, sep = "\t", header = T, as.is = T)
colnames(tx) = c("tchr", "tbeg", "tend", "tsrd", "qchr", "qbeg", "qend", "qsrd", "cid", "lev")
#tx = tx[tx$lev == 1,]

fh = file.path(ccfg$cdir, "../51_ortho/01.ortho.tbl")
th = read.table(fh, sep = "\t", header = T, as.is = T)[,1:8]
thsr = th[!is.na(th$cid) & (th$slen/th$tlen >= 0.5 | th$slen/th$qlen >= 0.5),]
gb = group_by(thsr, qid)
ths = summarise(gb, i = which(slen==max(slen))[1], idx = idx[i], tid = tid[i], tlen = tlen[i], qlen = qlen[i], slen = slen[i], cid = cid[i], lev = lev[i])
ths = ths[,c(3,4:5,1,6:9)]
ids_del_fus = thsr$tid[!thsr$idx %in% ths$idx]

fv = file.path(ccfg$cdir, "../31_sv/05.stb")
tv = read.table(fv, sep = "\t", header = T, as.is = T)[,c(2:5,7,8:11)]
tdel = tv[tv$tlen >= 50,]
tins = tv[tv$qlen >= 50,]

### tgt del
fg = file.path(tcfg$dir, "51.tbl")
tg = read.table(fg, sep = "\t", header = F, as.is = T)
colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
tg = tg[tg$type == 'cds',]
grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gb = group_by(tg, id)
dg = summarise(gb, len = sum(end-beg+1), chr = chr[1], beg = min(beg), end = max(end), fam = fam[1])
tdg = dg

gv = with(tdel, GRanges(seqnames = tchr, ranges = IRanges(tbeg+1, end = tend-1)))
tt = intersect_cds_sv(grg, gv, tg$id)
to = merge(dg, tt, by.x = "id", by.y = "idx")
to = cbind(to, pct_del = to$olen / to$len)
#tx = merge(th, to[,c('id', 'pct_del')], by.x = 'tid', by.y = 'id')

ids_del_1 = to$id[to$pct_del >= 0.5]
ids_del_2 = to$id[to$pct_del >= 0.2 & to$pct_del < 0.5 & !to$id %in% ths$tid]
ids_del = unique(c(ids_del_fus, ids_del_1, ids_del_2))

### qry ins
fg = file.path(qcfg$dir, "51.tbl")
tg = read.table(fg, sep = "\t", header = F, as.is = T)
colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
tg = tg[tg$type == 'cds',]
grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gb = group_by(tg, id)
dg = summarise(gb, len = sum(end-beg+1), chr = chr[1], beg = min(beg), end = max(end), fam = fam[1])
qdg = dg

gv = with(tins, GRanges(seqnames = qchr, ranges = IRanges(qbeg+1, end = qend-1)))
tt = intersect_cds_sv(grg, gv, tg$id)
to = merge(dg, tt, by.x = "id", by.y = "idx")
to = cbind(to, pct_ins = to$olen / to$len)

ids_ins_1 = to$id[to$pct_ins >= 0.5]
ids_ins_2 = to$id[to$pct_ins >= 0.2 & to$pct_ins < 0.5 & !to$id %in% ths$qid]
ids_ins = unique(c(ids_ins_1, ids_ins_2))

### unify
ty = ths
ty$qid[ty$tid %in% ids_del] = '-'
ty$tid[ty$qid %in% ids_ins] = "-"

tdgs = tdg[!tdg$id %in% unique(ty$tid),]
ttgt = data.frame(idx = NA, tid = tdgs$id, tlen = tdgs$len, qid = '', qlen = NA, slen = NA, cid = NA, lev = NA, stringsAsFactors = F)
tz = rbind(ty, ttgt)
tz$qid[tz$tid %in% ids_del] = '-'

qdgs = qdg[!qdg$id %in% unique(tz$qid),]
tqry = data.frame(idx = NA, tid = '', tlen = NA, qid = qdgs$id, qlen = qdgs$len, slen = NA, cid = NA, lev = NA, stringsAsFactors = F)
tz = rbind(tz, tqry)
tz$tid[tz$qid %in% ids_ins] = '-'

qst = tz$qid
qst[qst != '' & qst != '-'] = 'x'
tst = tz$tid
tst[tst != '' & tst != '-'] = 'x'

table(qst, tst)
stopifnot(nrow(tdg) == sum(tst=='x'))
stopifnot(nrow(qdg) == sum(qst=='x'))

tu = cbind(tz, tst = tst, qst = qst)
tu = tu[order(tu$idx), -1]
fo = file.path(ccfg$cdir, "../51_ortho/03.ortho.indel.tbl")
write.table(tu, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

### determine ins location
tj = to[to$id %in% ids_ins,]
tk = merge(tj, cbind(tins[,c('tchr','tbeg')], qidx = 1:nrow(tins)), by = 'qidx')

get_gene_up_pos <- function(rw, dg) {
  chr = rw['tchr']; pos = as.numeric(rw['tbeg'])
  dgs = dg[dg$chr == chr & dg$beg < pos & dg$end < pos,]
  if(nrow(dgs) == 0) {
    1
  } else {
    dgs = dgs[order(dgs$beg, decreasing = T),]
    dgs$end[1] + 1
  }
}
poss = apply(tk[,c('tchr','tbeg')], 1, get_gene_up_pos, tdg)

to = data.frame(chr = tk$tchr, pos = poss, tk[,c('id','len')])
fo = file.path(ccfg$cdir, "../51_ortho/04.ins.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')


cat(qname,"\n")
}

#### prepare to run mcl
tcfg = get_genome_cfg(tname)
ccfgs = get_comp_cfg(tname, qnames)
lst = list()
for (qname in qnames) {
  ccfg = ccfgs[[qname]]

  fi = file.path(ccfg$cdir, "../51_ortho/04.ins.tbl")
  ti = read.table(fi, header = T, sep = "\t", as.is = T)
  lst[[qname]] = cbind(ti, org = qname)
}
ds = do.call(rbind.data.frame, lst)

gb = group_by(ds, chr, pos)
to = summarise(gb, note = paste(sprintf("%s-%s-%d", org, id, len), collapse = " "))
to = to[order(to$chr, to$pos),]
to = cbind(idx = 1:nrow(to), to)

fo = file.path(Sys.getenv("misc3"), "comp.ortho.ins", "01.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

### run comp.ortho.ins.tbl to generate 05.group.tbl
##### combine single-ins with multi-ins
fi = file.path(dirw, "01.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
nids = as.numeric(sapply(ti$note, function(x) length(strsplit(x, split = " ")[[1]])))
ti = ti[nids == 1, ]
tsg = within(ti, {
  qry = sapply(strsplit(note, "-"), function(x) x[[1]][1])
  qid = sapply(strsplit(note, "-"), function(x) x[[2]][1])
  rm(note)
})

tw2 = reshape(tsg[,c(-2,-3)], direction = 'wide', timevar = 'qry', 
  idvar = 'idx', times = list(qnames = qnames))
tw2 = tw2[,c('idx', paste("qid",qnames,sep="."))]
tw2 = merge(tsg[,1:3], tw2, by = 'idx')
tw2 = cbind(tw2[,1:3], tid = '-', tw2[4:ncol(tw2)])


fi = file.path(Sys.getenv("misc3"), "comp.ortho.ins", "05.group.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti = ti[,c(colnames(ti)[1:3], qnames)]
colnames(ti)[4:ncol(ti)] = paste("qid", colnames(ti)[4:ncol(ti)], sep = ".")
tw3 = cbind(ti[,1:3], tid = '-', ti[4:ncol(ti)])
#tw3$idx = 1:nrow(tw3)

stopifnot(colnames(tw2) == colnames(tw3))
tw = rbind(tw2, tw3)
tw = tw[order(tw$idx),]
tw[tw == ''] = NA
colnames(tw)[4:ncol(tw)] = c(tname, qnames)

##### fill missing data (use 'tw' from above step)
fg = file.path(tcfg$dir, "51.tbl")
tg = read.table(fg, sep = "\t", header = F, as.is = T)
colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
tg = tg[tg$type == 'cds',]
grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gb = group_by(tg, id)
dg = summarise(gb, len = sum(end-beg+1), chr = chr[1], beg = min(beg), end = max(end), fam = fam[1])

dz = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
colnames(dz) = c("chr", "size")

fi = file.path(dirw, "01.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)

get_gene_down_pos <- function(rw, dg, dz) {
  chr = rw['chr']; pos = as.numeric(rw['pos'])
  dgs = dg[dg$chr == chr & dg$beg > pos & dg$end > pos,]
  if(nrow(dgs) == 0) {
    dz$size[dz$chr == chr]
  } else {
    dgs = dgs[order(dgs$beg),]
    dgs$end[1]
  }
}
poss = apply(ti, 1, get_gene_down_pos, dg, dz)
ti = cbind(ti[,1:3], end = poss)
gri = with(ti, GRanges(seqnames = chr, ranges = IRanges(pos, end = end)))


intersect_fill_ins <- function(gr1, gr2) {
  t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
  t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), score = mcols(gr2)$score, stringsAsFactors = F)
  
  fbd1 = 'xtest1.bed'
  fbd2 = 'xtest2.bed'
  fres = 'xout.bed'
  options(scipen = 999)
  write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
  options(scipen = 0)
  system(sprintf("intersectBed -wao -a %s -b %s > %s", fbd1, fbd2, fres))

  t3 = read.table(fres, sep = "\t", header = F, as.is = T)
  colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'cid', 'olen')
  system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

  gp = group_by(t3, idx)
  t4 = dplyr::summarise(gp, ncid = length(unique(cid)), obeg = min(qbeg), oend = max(qend))
  t5 = merge(t1, t4, by = 'idx')
  stopifnot(nrow(t1) == nrow(t5))
  apply(t5, 1, function(x) as.numeric(x['ncid']) == 1 & as.numeric(x['obeg']) <= as.numeric(x['beg']) & as.numeric(x['oend']) >= as.numeric(x['end']))
}

qname = qnames[1]
for (qname in qnames) {
  ccfg = ccfgs[[qname]]
  fx = file.path(ccfg$cdir, "31.9/gal")
  tx = read.table(fx, sep = "\t", header = T, as.is = T)[,1:6]
  grx = with(tx, GRanges(seqnames = tId, ranges = IRanges(tBeg, end = tEnd), score = id))
  status = intersect_fill_ins(gri, grx)
  tws = merge(tw[,c('idx', qname)], cbind(ti, sta = status), by = 'idx')
  sta = tws[,qname]
  sta[is.na(sta) & tws$sta] = '-'
  tw[,qname] = sta
}

fo = file.path(dirw, "08.tbl")
write.table(tw, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')



