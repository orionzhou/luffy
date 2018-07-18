require(plyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(reshape2)
require(RColorBrewer)
require(ape)
require(bios2mds)
require(gridBase)
require(colorRamps)
source("comp.fun.R")

dirw = file.path(Sys.getenv("misc3"), "comp.ortho")

### run comp.ortho.ins.R to generate comp.ortho.ins/08.tbl
#### create 01.tbl   HM101+del+ins ortho-map
fg = file.path(tcfg$dir, "51.tbl")
tg = read.table(fg, sep = "\t", header = F, as.is = T)
colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
tg = tg[tg$type == 'mrna', c('chr','beg','id')]

tl = data.frame()
for (qname in qnames) {
  ccfg = ccfgs[[qname]]

  fi = file.path(ccfg$cdir, "../51_ortho/03.ortho.indel.tbl")
  ti = read.table(fi, header = T, sep = "\t", as.is = T)
  tls = cbind(qry = qname, ti[ti$tst == 'x', c('tid','qid')])
  tl = rbind(tl, tls)
}
tw1 = reshape(tl, direction = 'wide', timevar = 'qry', idvar = 'tid',
  times = list(qnames = qnames))
tw1 = tw1[,c('tid', paste("qid", qnames, sep="."))]
tw1 = merge(tg, tw1, by.x = 'id', by.y = 'tid')
colnames(tw1)[1:3] = c('tid','chr','pos')
tw1 = cbind(idx = NA, tw1[,c(2:3,1,4:ncol(tw1))])
colnames(tw1)[4:ncol(tw1)] = c(tname, qnames)


fi = file.path(Sys.getenv("misc3"), "comp.ortho.ins", "08.tbl")
tw2 = read.table(fi, header = T, sep = "\t", as.is = T)

stopifnot(colnames(tw1) == colnames(tw2))
tw = rbind(tw1, tw2)
tw = tw[order(tw$chr, tw$pos, tw$idx),]
tw$idx = 1:nrow(tw)
norgs = apply(tw[,orgs], 1, function(x) sum(!x %in% c('','-')))
table(norgs)

fo = file.path(dirw, "01.tbl")
write.table(tw, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')


##### create 05.ortho.tbl and 06.no.ortho.tbl 
fi = file.path(dirw, "01.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
norgs = apply(ti[,orgs], 1, function(x) sum(!x %in% c('','-')))
table(norgs)

cnames = colnames(ti)
lst = apply(ti, 1, function(x, cnames) {
  idx = which(! x[4:length(cnames)] %in% c('', '-'))[1] + 4 - 1
  list('org'=cnames[idx], 'gid'=x[idx])
}, cnames)
to = as.data.frame(do.call(rbind, lst), stringsAsFactors = F)
to = data.frame(org = as.character(to[,1]), gid = as.character(to[,2]), idx = ti$idx)

fo = file.path(dirw, "05.ortho.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')


til = melt(ti, id.vars = c('idx','chr','pos'), measure.vars = orgs,
  variable.name = 'org', value.name = 'gid')
tort = til[til$gid != '' & til$gid != '-',]

tg = data.frame()
for (org in orgs) {
  cfg = cfgs[[org]]
  fg = file.path(cfg$dir, "51.gtb")
  gids = read.table(fg, header = T, sep = "\t", as.is = T)[,1]
  tgs = data.frame(org = org, gid = gids, stringsAsFactors = F)
  tg = rbind(tg, tgs)
}
tm = merge(tg, tort[,c('org','gid','idx')], by = c("org", "gid"), all = T)
stopifnot(nrow(tm) == nrow(tg))

tsg = tm[is.na(tm$idx),]
fo = file.path(dirw, "06.no.ortho.tbl")
write.table(tsg[,1:2], fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')


### run comp.ortho.pl to create 09.tbl (blat no.ortho.fas against ortho.fas)
##### fill in 01.tbl with 09.tbl blat result
fw = file.path(dirw, "01.tbl")
tw = read.table(fw, header = T, sep = "\t", as.is = T)
stopifnot(sum(tw$idx - 1:nrow(tw)) == 0)

tl = melt(tw, id.vars = c('idx','chr','pos'), measure.vars = orgs,
  variable.name = 'org', value.name = 'gid')
tsg = tl[tl$gid == '',]


fi = file.path(dirw, "09.tbl")
ti = read.table(fi, sep = "\t", header = F, stringsAsFactors = F)
colnames(ti) = c("tid", "qid", "tcov", "qcov", "score")
ti = within(ti, {
  idx = as.numeric(sapply(strsplit(tid, "-"), function(x) x[[3]][1]))
  org = sapply(strsplit(qid, "-"), function(x) x[[1]][1])
  gid = sapply(strsplit(qid, "-"), function(x) x[[2]][1])
  rm(tid, qid)
})

tid = tw[,c('idx',orgs)]
tlo = tw[,1:3]
tst = tid[,-1]
tst[tst != '' & tst != '-'] <- 'syn'
tst = cbind(idx = tid$idx, tst)

qname = 'HM340'
for (qname in qnames) {
  tss = tl[tl$gid == '' & tl$org == qname,]
  tis = ti[ti$idx %in% tss$idx & ti$org == qname,]
  tis = tis[order(tis$idx),]
  
  gids = unique(tis$gid)
  hs = rep(0, length(gids))
  names(hs) = gids
  
  x = Rle(tis$idx)
  ibeg = start(x); iend = end(x); idxs = runValue(x)
  ogids = c(); oidxs = c()
  for (i in 1:length(ibeg)) {
    ds = tis[ibeg[i]:iend[i],]
    ds = cbind(ds, sta = as.numeric(hs[ds$gid]))
    ds = ds[ds$sta == 0,]
    if(nrow(ds) > 0) {
      ds = ds[order(ds$score, decreasing = T),]
      hs[ds$gid[1]] = 1
      ogid = ds$gid[1]
      
      oidxs = c(oidxs, idxs[i])
      ogids = c(ogids, ogid)
    }
  }
  stopifnot(sum(tid[oidxs, qname] != '') == 0)
  tid[oidxs, qname] = ogids
  tst[oidxs, qname] = 'rbh'
  
  
  tss = tl[tl$gid == '-' & tl$org == qname,]
  tis = ti[ti$idx %in% tss$idx & ti$org == qname,]
  tis = tis[order(tis$idx),]
  
  gids = unique(tis$gid)
  hs = rep(0, length(gids))
  names(hs) = gids
  
  x = Rle(tis$idx)
  ibeg = start(x); iend = end(x); idxs = runValue(x)
  ogids = c(); oidxs = c()
  for (i in 1:length(ibeg)) {
    ds = tis[ibeg[i]:iend[i],]
    ds = cbind(ds, sta = as.numeric(hs[ds$gid]))
    ds = ds[ds$sta == 0,]
    if(nrow(ds) > 0) {
      ds = ds[order(ds$score, decreasing = T),]
      hs[ds$gid[1]] = 1
      ogid = ds$gid[1]
      
      oidxs = c(oidxs, idxs[i])
      ogids = c(ogids, ogid)
    }
  }
  stopifnot(sum(tid[oidxs, qname] != '-') == 0)
  tid[oidxs, qname] = ogids
  tst[oidxs, qname] = '-|rbh'
  
  
  x = table(tst[,qname])
  cat(qname, x['syn'], x['rbh'], x['-|rbh'], x['-'], x[1], "\n")
}
stopifnot(sum(tid=='')==sum(tst==''))
norgs = apply(tst[,orgs], 1, function(x) sum(!x %in% c('','-')))
table(norgs)

fo = file.path(dirw, "11.loc.tbl")
write.table(tlo, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')
fo = file.path(dirw, "11.gid.tbl")
write.table(tid, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')
fo = file.path(dirw, "11.sta.tbl")
write.table(tst, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

## output no-ortho
fi = file.path(dirw, "11.gid.tbl")
tid = read.table(fi, header = T, sep = "\t", as.is = T)

tg = data.frame()
for (org in orgs) {
  cfg = cfgs[[org]]
  fg = file.path(cfg$dir, "51.gtb")
  gids = read.table(fg, header = T, sep = "\t", as.is = T)[,1]
  sgids = gids[!gids %in% unique(tid[,org])]
  if(org == tname) {
    stopifnot(length(sgids) == 0)
  } else {
    tgs = data.frame(org = org, gid = sgids, stringsAsFactors = F)
    tg = rbind(tg, tgs)
    cat(org, nrow(tgs), "\n")
  }
}

fo = file.path(dirw, "15.no.ortho.tbl")
write.table(tg, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

## run comp.ortho.pl to generate 17.group.tbl
##### add grouped-RBH ortho
fi = file.path(dirw, "11.gid.tbl")
tid = read.table(fi, header = T, sep = "\t", as.is = T)
fi = file.path(dirw, "11.loc.tbl")
tlo = read.table(fi, header = T, sep = "\t", as.is = T)
fi = file.path(dirw, "11.sta.tbl")
tst = read.table(fi, header = T, sep = "\t", as.is = T)
stopifnot(sum(tid$idx - 1:nrow(tid)) == 0)
stopifnot(sum(tlo$idx - 1:nrow(tid)) == 0)
stopifnot(sum(tst$idx - 1:nrow(tid)) == 0)

fi = file.path(dirw, "17.group.tbl")
tgr = read.table(fi, header = T, sep = "\t", as.is = T)
stopifnot(sum(!is.na(tgr$HM101)) == 0)
tgr$HM101 = ''
tgr = tgr[,orgs]

norgs = apply(tgr, 1, function(x) sum(x[2:length(x)] != ''))
table(norgs)

fi = file.path(dirw, "15.no.ortho.tbl")
tn = read.table(fi, header = T, sep = "\t", as.is = T)
tsg = data.frame()
for (org in qnames) {
  ids_no = tn$gid[tn$org == org]
  ids_gr = unique(tgr[,org])
  ids_sg = ids_no[!ids_no %in% ids_gr]
  cat(org, length(ids_sg), "\n")
  tz = matrix("", ncol = length(orgs), nrow = length(ids_sg))
  colnames(tz) = orgs
  tz[,org] = ids_sg
  tsg = rbind(tsg, data.frame(tz, stringsAsFactors = F))
}
tz = rbind(tgr, tsg)

ctid = cbind(idx = nrow(tid)+(1:nrow(tz)), tz)
ctst = ctid[,2:ncol(ctid)]
ctst[ctst != ''] = 'rbh'
ctst = cbind(idx = ctid$idx, ctst)
ctlo = data.frame(idx = nrow(tid)+(1:nrow(tz)), chr = 'chrZ', pos = 10000*(1:nrow(ctid)), stringsAsFactors = F)

tid = rbind(tid, ctid); tlo = rbind(tlo, ctlo); tst = rbind(tst, ctst)
fo = file.path(dirw, "21.loc.tbl")
write.table(tlo, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')
fo = file.path(dirw, "21.gid.tbl")
write.table(tid, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')
fo = file.path(dirw, "21.sta.tbl")
write.table(tst, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')
norgs = apply(tst[,orgs], 1, function(x) sum(!x %in% c('','-')))
table(norgs)


tl = melt(tsg, id.vars = NULL, measure.vars = orgs,
  variable.name = 'org', value.name = 'gid')
tl = tl[tl$gid != '',]

fo = file.path(dirw, "29.singleton.tbl")
write.table(tl, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

### run comp.ortho.pl to create 22.cat.tbl

##### generate final score matrix
fi = file.path(dirw, "21.gid.tbl")
tid = read.table(fi, header = T, sep = "\t", as.is = T)
fi = file.path(dirw, "21.loc.tbl")
tlo = read.table(fi, header = T, sep = "\t", as.is = T)
fi = file.path(dirw, "21.sta.tbl")
tst = read.table(fi, header = T, sep = "\t", as.is = T)
stopifnot(sum(tid$idx - 1:nrow(tid)) == 0)
stopifnot(sum(tlo$idx - 1:nrow(tid)) == 0)
stopifnot(sum(tst$idx - 1:nrow(tid)) == 0)
norgs = apply(tst, 1, function(x) sum(x[2:length(x)] %in% c('syn','rbh')))
table(norgs)

get_pw_comps <- function(ids) {
  comps = c()
  for (i in 1:(length(ids)-1)) {
    id1 = ids[i]
    for (j in (i+1):length(ids)) {
      id2 = ids[j]
      comps = c(comps, sprintf("%s-%s", id1, id2))
    }
  }
  comps
}
comps = get_pw_comps(orgs)

get_pd <- function(idx, norgs, orgs, comps, dira, get_pw_comps) {
  res = rep(NA, length(comps))
  names(res) = comps
  norg = norgs[idx]
  if(norg > 1) {
    fa = sprintf("%s/%d.fas", dira, idx)
    aln = import.fasta(fa)
    dis = mat.dis(aln, aln)
    orgs1 = sapply(strsplit(rownames(dis), "[-]"), "[[", 1)
    rownames(dis) = orgs1; colnames(dis) = orgs1
    corgs = orgs[orgs %in% orgs1]
    dis = dis[corgs, corgs]
    ccomps = get_pw_comps(corgs)
    cvals = dis[lower.tri(dis)]
    res[ccomps] = cvals
  }
  res
}
dira = file.path(dirw, "25_aln")

cl = makeCluster(detectCores())
cluster_fun <- function() { 
  require(bios2mds)
}
clusterCall(cl, cluster_fun)

ptm <- proc.time()
#y = sapply(1:5, get_pd, norgs, orgs, comps, dira, get_pw_comps)
#y = parSapply(cl, 1:1000, get_pd, norgs, orgs, comps, dira, get_pw_comps)
y = parSapply(cl, 1:nrow(tid), get_pd, norgs, orgs, comps, dira, get_pw_comps)
proc.time() - ptm

to = t(y)
fo = file.path(dirw, "28.dist.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')
