require(GenomicRanges)
require(dplyr)
options(scipen=999)

locCluster <- function(pos, wsize) {
  npos = length(pos)
  df = data.frame(id=1:npos, pos=pos, cluster=1:npos)
  df = df[order(df$pos),]
  for (i in 1:npos) {
    for (j in (i+1):npos) {
      if(j>npos) {
        next
      }
      if(df$pos[j] - df$pos[i] <= wsize) {
        df$cluster[j] = df$cluster[i]
      } else {
        next
      }
    }
  }
  clusters = unique(df$cluster)
  tmp = cbind(cluster1=clusters, cluster2=1:length(clusters))
  x = merge(df, tmp, by.x='cluster', by.y='cluster1')
  
  df2 = data.frame(id=x$id, pos=x$pos, cluster=x$cluster2, cluster_y=1)
  df2 = df2[order(df2$pos),]
  
  clusterP = ''
  for (i in 1:npos) {
    if(df2$cluster[i] == clusterP) {
      df2$cluster_y[i] = df2$cluster_y[i-1] + 1
    } else {
      clusterP = df2$cluster[i]
    }
  }

#  hist(as.matrix(table(df2$cluster)), xlab='cluster size', main=paste(wsize, 'bp', sep=''))
  df2
}

locStr2List <- function(str) {
  tmp = gregexpr("complement([[:graph:]]+)", str)[[1]]
  strand = 1
  if(tmp[1] == 1) { strand = -1 }
  tmp = gregexpr("([0-9]+\\.\\.[0-9]+)", str)[[1]]
  strs = substring(str, tmp, tmp+attr(tmp, "match.length")-1)
  tmp = strsplit(strs, "\\.\\.")
  begs = as.numeric(sapply(tmp, "[", 1))
  ends = as.numeric(sapply(tmp, "[", 2))
  list(strand=strand, begs=begs, ends=ends)
}

locStr2Df <- function(str, seqid='chrN', type='') {
  tmp = gregexpr("complement([[:graph:]]+)", str)[[1]]
  strand = 1 
  if(tmp[1] == 1) { strand = -1 }
  tmp = gregexpr("([0-9]+\\.\\.[0-9]+)", str)[[1]]
  if(tmp[1] == -1) {
    data.frame()
  } else {
    strs = substring(str, tmp, tmp+attr(tmp, "match.length")-1)
    tmp = strsplit(strs, "\\.\\.")
    begs = as.numeric(sapply(tmp, "[", 1))
    ends = as.numeric(sapply(tmp, "[", 2))
    data.frame(chr=seqid, beg=begs, end=ends, type=type)
  }
}
get_loc_gene <- function(g) {
  id = as.character(g['id'])
  chr = as.character(g['chr'])
  locC = as.character(g['locC'])
  locI = as.character(g['locI'])
  loc5 = as.character(g['loc5'])
  loc3 = as.character(g['loc3'])
  df = rbind(locStr2Df(locC,chr,"cds"), locStr2Df(locI,chr,'intron'), locStr2Df(loc5,chr,'utr5'), locStr2Df(loc3,chr,'utr3'))
  cbind(df, id=id)
}

intersect_basepair <- function(gr1, gr2) {
  t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
  gr2 = reduce(gr2)
  t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), stringsAsFactors = F)
  
  fbd1 = 'xtest1.bed'
  fbd2 = 'xtest2.bed'
  fres = 'xout.bed'
  options(scipen = 999)
  write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
  options(scipen = 0)
  system(sprintf("intersectBed -wao -a %s -b %s > %s", fbd1, fbd2, fres))

  t3 = read.table(fres, sep = "\t", header = F, as.is = T)
  colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
  system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

  gp = group_by(t3, idx)
  t4 = dplyr::summarise(gp, olen = sum(olen))
  as.numeric(t4$olen)
}
intersect_score <- function(gr1, gr2) {
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
  colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'score', 'olen')
  t3$score[t3$olen == 0] = NA
  system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

  gp = group_by(t3, idx)
  t4 = dplyr::summarise(gp, score = sum(score, na.rm = T))
  as.numeric(t4$score)
}
intersect_count <- function(gr1, gr2) {
  t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
  t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), stringsAsFactors = F)
  
  fbd1 = 'xtest1.bed'
  fbd2 = 'xtest2.bed'
  fres = 'xout.bed'
  options(scipen = 999)
  write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
  options(scipen = 0)
  system(sprintf("intersectBed -wao -a %s -b %s > %s", fbd1, fbd2, fres))

  t3 = read.table(fres, sep = "\t", header = F, as.is = T)
  colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
  system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

  gp = group_by(t3, idx)
  t4 = dplyr::summarise(gp, cnt = sum(qidx != '.'))
  as.numeric(t4$cnt)
}
intersect_idx <- function(gr1, gr2) {
  t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
  t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), stringsAsFactors = F)
  
  fbd1 = 'xtest1.bed'
  fbd2 = 'xtest2.bed'
  fres = 'xout.bed'
  options(scipen = 999)
  write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
  options(scipen = 0)
  system(sprintf("intersectBed -wao -a %s -b %s > %s", fbd1, fbd2, fres))

  t3 = read.table(fres, sep = "\t", header = F, as.is = T)
  colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
  system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

  t3[t3$qidx != '.', c('chr','beg','end','idx','qidx','olen')]
}
intersect_idx_maxovlp <- function(gr1, gr2) {
  t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
  t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), stringsAsFactors = F)
  
  fbd1 = 'xtest1.bed'
  fbd2 = 'xtest2.bed'
  fres = 'xout.bed'
  options(scipen = 999)
  write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
  options(scipen = 0)
  system(sprintf("intersectBed -wao -a %s -b %s > %s", fbd1, fbd2, fres))

  t3 = read.table(fres, sep = "\t", header = F, as.is = T)
  colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
  system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

  gp = group_by(t3, idx)
  t4 = dplyr::summarise(gp, qidx = qidx[which(olen == max(olen))[1]])
  as.numeric(t4$qidx)
}
intersect_idx_t <- function(gr1, gr2) {
  t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
  t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), stringsAsFactors = F)
  
  fbd1 = 'xtest1.bed'
  fbd2 = 'xtest2.bed'
  fres = 'xout.bed'
  options(scipen = 999)
  write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
  options(scipen = 0)
  system(sprintf("intersectBed -wo -a %s -b %s > %s", fbd1, fbd2, fres))

  t3 = read.table(fres, sep = "\t", header = F, as.is = T)
  colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
  system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

  as.numeric(t3$idx)
}
intersect_cds_sv <- function(gr1, gr2, cdsidx) {
  t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = cdsidx, stringsAsFactors = F)
  t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), stringsAsFactors = F)
  
  fbd1 = 'xtest1.bed'
  fbd2 = 'xtest2.bed'
  fres = 'xout.bed'
  options(scipen = 999)
  write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
  options(scipen = 0)
  system(sprintf("intersectBed -wao -a %s -b %s > %s", fbd1, fbd2, fres))

  t3 = read.table(fres, sep = "\t", header = F, as.is = T)
  colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
  system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

  gp = group_by(t3, idx, qidx)
  t4 = dplyr::summarise(gp, olen = sum(olen))
  gp2 = group_by(t4, idx)
  t5 = dplyr::summarise(gp2, qidx = qidx[which(olen == max(olen))][1], olen = sum(olen))
  t5
}
intersect_gene <- function(gr1, gr2, cdsidx) {
  t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
  t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = cdsidx, stringsAsFactors = F)
  
  fbd1 = 'xtest1.bed'
  fbd2 = 'xtest2.bed'
  fres = 'xout.bed'
  options(scipen = 999)
  write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
  options(scipen = 0)
  system(sprintf("intersectBed -wo -a %s -b %s > %s", fbd1, fbd2, fres))

  t3 = read.table(fres, sep = "\t", header = F, as.is = T, quote = "")
  colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
  system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

  gp = group_by(t3, idx, qidx)
  t4 = dplyr::summarise(gp, olen = sum(olen))
  t4
}

