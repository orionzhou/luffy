require(dplyr)
require(plyr)
require(rtracklayer)
#require(rbamtools)
require(Rsamtools)
require(seqminer)

read_seqinfo <- function(fsize) {
  tsize = read.table(fsize, header = F, sep = "\t",  as.is = T, 
    col.names = c("id", "size"))
  Seqinfo(tsize$id, seqlengths = tsize$size)
}
get_genome_cfg <- function(org, tname = "HM101") {
  org = toupper(org)
  gdir = file.path('/home/youngn/zhoux379/data/genome', org)
  
  size = file.path(gdir, "15.sizes")
  seqinfo = read_seqinfo(file.path(gdir, "15.sizes"))
  
  gap = file.path(gdir, "16.gap.bed")
  gapz = file.path(gdir, "16.gap.bed.gz")
  
  gene = file.path(gdir, "51.tbl")
  genez = file.path(gdir, "51.tbl.gz")

  dna = file.path(gdir, '11_genome.fas')
  protein = file.path(gdir, '51.fas')
  
  mapp = file.path(gdir, "18_stat_k60/15_mapp.bw")

  cfg = list(gdir = gdir, size = size, seqinfo = seqinfo, 
    gap = gap, gapz = gapz, 
    gene = gene, genez = genez,
    dna = dna, protein = protein, mapp = mapp)
  
  orgp = gsub(".AC", "", org)
  f_pacbio = sprintf("%s/pacbio/%s_%s/15.bam", Sys.getenv("misc3"), orgp, org)
  #if(file.exists(f_pacbio)) cfg[['pacbio']] = bamReader(f_pacbio, idx = T)

  f_rnaseq = sprintf("%s/rnaseq/mt/22_tophat/%s_%s/accepted_hits.bam", Sys.getenv("misc2"), orgp, org)
  #if(file.exists(f_rnaseq)) cfg[['rnaseq']] = bamReader(f_rnaseq, idx = T)
  
  if(org != tname) {
    cdir = sprintf("%s/%s_%s/23_blat", '/home/youngn/zhoux379/data/misc3', org, tname)
    tgal = sprintf("%s/31.9/gal.gz", cdir)
    qgal = sprintf("%s/41.9/gal.gz", cdir)
    tgax = sprintf("%s/31.9/gax.gz", cdir)
    qgax = sprintf("%s/41.9/gax.gz", cdir)
    tsnp = sprintf("%s/31.9/snp.gz", cdir)
    qsnp = sprintf("%s/41.9/snp.gz", cdir)
    fpct = sprintf("%s/31.9/pct.bw", cdir)
    
    vdir = sprintf("%s/hapmap_mt40/12_ncgr", '/home/youngn/zhoux379/data/misc3')
    fsnp = sprintf("%s/44_snp/%s.tbl.gz", vdir, org)
    fcov = sprintf("%s/35_cov/%s.bw", vdir, org)
    fcovab = sprintf("%s/36_abcov/%s.bw", vdir, org)
    
    ccfg = list(cdir = cdir, 
      tgal = tgal, tgax = tgax, tsnp = tsnp, qsnp = qsnp, tpct = fpct,
      vdir = vdir, vsnp = fsnp, vcov = fcov)
    cfg = c(cfg, ccfg)
  }
  cfg
}
get_genome_cfgs <- function(orgs, tname = "HM101") {
  cfgs = list()
  for (org in orgs) {
    cfgs[[org]] = get_genome_cfg(org)
  }
  cfgs
}
granges2df <- function(gr) {
  ds = data.frame(chr = as.character(seqnames(gr)), beg = start(gr), 
    end = end(gr), srd = as.character(strand(gr)), stringsAsFactors = F)
  ds$srd[ds$srd == "*"] = "+"
  ds
}
get_genome_composition <- function(org, utr.merge = F) {
  dirw = file.path('/home/youngn/zhoux379/data/genome', org)
  fsize = file.path(dirw, "15.sizes")
  fgap = file.path(dirw, "16.gap.bed")
  stopifnot(file.exists(fsize) & file.exists(fgap))
  
  tlen = read.table(fsize, sep = "\t", header = F, as.is = T)
  grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))

  tgap = read.table(fgap, sep = "\t", header = F, as.is = T)
  grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

  gra = GenomicRanges::setdiff(grt, grp)

  f1 = file.path(dirw, "51.tbl")
  stopifnot(file.exists(f1))
  
  t1 = read.table(f1, header = F, sep = "\t", as.is = T)
  colnames(t1) = c("chr", "beg", "end", "srd", "id", "type", "fam")
  
  tt = t1[t1$type == 'cds',]
  gcds = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
  gcds = GenomicRanges::intersect(gra, reduce(gcds))

  tt = t1[t1$type == 'intron',]
  gito = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
  gito = GenomicRanges::intersect(GenomicRanges::setdiff(gra, gcds), reduce(gito))

  tt = t1[t1$type == 'utr5',]
  gut5 = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
  gut5 = GenomicRanges::intersect(GenomicRanges::setdiff(gra, c(gcds, gito)), reduce(gut5))

  tt = t1[t1$type == 'utr3',]
  gut3 = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
  gut3 = GenomicRanges::intersect(GenomicRanges::setdiff(gra, c(gcds, gito,   gut5)), reduce(gut3))

  gutr = GenomicRanges::union(gut5, gut3)
  gitr = GenomicRanges::setdiff(gra, c(gcds, gito, gut5, gut3))

  if(utr.merge) {
    rbind( cbind(granges2df(gcds), type = 'cds'),
      cbind(granges2df(gito), type = 'intron'),
      cbind(granges2df(gutr), type = 'utr'),
      cbind(granges2df(gitr), type = 'intergenic')
    )
  } else {
    rbind( cbind(granges2df(gcds), type = 'cds'),
      cbind(granges2df(gito), type = 'intron'),
      cbind(granges2df(gut5), type = 'utr5'),
      cbind(granges2df(gut3), type = 'utr3'),
      cbind(granges2df(gitr), type = 'intergenic')
    )
  }
}

trim_gax <- function(ds, beg, end) {
  for (i in 1:nrow(ds)) {
    if(ds$tbeg[i] < beg) {
      if(as.character(ds$tsrd[i]) == as.character(ds$qsrd[i])) {
        ds$qbeg[i] = ds$qbeg[i] + (beg - ds$tbeg[i])
      } else {
        ds$qend[i] = ds$qend[i] - (beg - ds$tbeg[i])
      }
      ds$tbeg[i] = beg
    }
    if(end < ds$tend[i]) {
      if(as.character(ds$tsrd[i]) == as.character(ds$qsrd[i])) {
        ds$qend[i] = ds$qend[i] - (ds$tend[i] - end)
      } else {
        ds$qbeg[i] = ds$qbeg[i] + (ds$tend[i] - end)
      }
      ds$tend[i] = end
    }
  }
  ds
}
read_gax_simple <- function(fgax, gr) {
  gr = reduce(gr)
  
  tg = data.frame()
  
  gax = open(TabixFile(fgax))
  grs = gr[as.character(seqnames(gr)) %in% seqnamesTabix(gax)]
  if(length(grs) == 0) return(NULL)
  x = scanTabix(gax, param = grs)
  close(gax)
  
  txts = rapply(x, c)
  if(length(txts) == 0) return(NULL)
  
  lc = list()
  for (i in 1:length(grs)) {
    if(length(x[[i]]) == 0) next
    df1 = parse_tabix(x[[i]])
    colnames(df1) = c('tid', 'tbeg', 'tend', 'tsrd', 
      'qid', 'qbeg', 'qend', 'qsrd', 'cid', 'lev')
    tgs = trim_gax(df1, start(grs)[i], end(grs)[i])
    tgs = cbind(tgs, len = tgs$qend - tgs$qbeg + 1)
    
    tg = rbind(tg, tgs)
  }
  tg
}
read_gal <- function(fgal, gr, minp = 0.01) {
  gr = reduce(gr)
  tc = data.frame()
  for (i in 1:length(gr)) {
  	locstr = sprintf("%s:%d-%d", seqnames(gr)[i], start(gr)[i], end(gr)[i])
    df1 = tabix.read.table(fgal, locstr)[,1:19]
    if(nrow(df1) == 0) next
  	colnames(df1) = c('cid', 'tid', 'tbeg', 'tend', 'tsrd', 'tsize',
      'qid', 'qbeg', 'qend', 'qsrd', 'qsize', 
      'lev', 'ali', 'mat', 'mis', 'qN', 'tN', 'ident', 'score')
    tcs = trim_gax(df1, start(gr)[i], end(gr)[i])
    tc = rbind(tc, tcs)
	}
  tc = tc[order(tc$tid, tc$tbeg, tc$tend),]
  
  idxs = tc$ali >= sum(tc$ali) * minp
  cids = tc$cid[idxs]
  list(tc = tc[idxs,])
}
read_gax <- function(fgax, gr, minp = 0.02) {
  gr = reduce(gr)
  tg = data.frame()
  tc = data.frame()
  
  lc = list()
  for (i in 1:length(gr)) {
  	locstr = sprintf("%s:%d-%d", seqnames(gr)[i], start(gr)[i], end(gr)[i])
    df1 = tabix.read.table(fgax, locstr)
    if(nrow(df1) == 0) next
    colnames(df1) = c('tid', 'tbeg', 'tend', 'tsrd', 
      'qid', 'qbeg', 'qend', 'qsrd', 'cid', 'lev')
    df1$cid = as.integer(df1$cid)
    tgs = trim_gax(df1, start(gr)[i], end(gr)[i])
    tgs = cbind(tgs, len = tgs$qend - tgs$qbeg + 1)
    
    gb = dplyr::group_by(tgs, cid)
    tcs = dplyr::summarise(gb, 
      tid = unique(tid), tbeg = min(tbeg), tend = max(tend), 
      tsrd = unique(tsrd),
      qid = unique(qid), qbeg = min(qbeg), qend = max(qend), 
      qsrd = unique(qsrd),
      ali = sum(len), gapo = length(len) - 1)
    for (cid in tcs$cid) {
      if(is.null(lc[[as.character(cid)]])) {
        lc[[as.character(cid)]] = 0
        next
      }
      lc[[as.character(cid)]] = lc[[as.character(cid)]] + 1
      ncid = cid + lc[[as.character(cid)]] / 100
      tcs$cid[tcs$cid == cid] = ncid
      tgs$cid[tgs$cid == cid] = ncid
    }
    tg = rbind(tg, tgs)
    tc = rbind(tc, tcs)
  }

  tc = cbind(tc, mis = 0)
  qgap = tc$qend - tc$qbeg + 1 - tc$ali
  tgap = tc$tend - tc$tbeg + 1 - tc$ali
  score = tc$ali * 1 + tc$mis * (-2) + tc$gapo * (-5) + 
    (qgap + tgap - tc$gapo) * (-2)
  tc = cbind(tc, score = score)
  tc = tc[order(tc$tid, tc$tbeg, tc$tend),]
  
  idxs = tc$ali >= sum(tc$ali) * minp
  cids = tc$cid[idxs]
  list(tg = tg[tg$cid %in% cids,], tc = tc[idxs,])
}

read_tabix <- function(ftbx, gr) {
  gr = reduce(gr)
  to = data.frame()
  for (i in 1:length(gr)) {
  	locstr = sprintf("%s:%d-%d", seqnames(gr)[i], start(gr)[i], end(gr)[i])
    df1 = tabix.read.table(ftbx, locstr)
    if(nrow(df1) == 0) next
    to = rbind(to, df1)
	}
  to
}
pileupReads <- function(dr){
    dr <- dr[order(dr$beg),]
    minbeg <- min(dr$beg); maxend <- max(dr$end)

    yread <- c(); ypos <- c()
    yread[1] <- minbeg - 1
    for (i in 1:nrow(dr)){
        beg <- dr$beg[i]
        placed <- F
        y <- 1
        while(!placed){
            if(yread[y] < beg){
                ypos[i] <- y
                yread[y] <- dr$end[i]
                placed <- T
            }
            y <- y + 1
            if(y > length(yread)){
                yread[y] <- minbeg - 1
            }
        }
    }
    ypos
}
read_bam <- function(bam, gr, pileup = F) {
  seqmap = getRefData(bam)
  
  res = data.frame()
  for (i in 1:length(gr)) {
    chr = as.character(seqnames(gr))[i]; beg = start(gr)[i]; end = end(gr)[i]
    refid = seqmap$ID[seqmap$SN == chr]
    reads = bamRange(bam, c(refid, beg, end))
    td <- as.data.frame(reads)
    if(is.null(td) | nrow(td) == 0) next
    tr = data.frame(chr = chr, beg = td$position + 1, end = NA, 
      len = nchar(td$seq), srd = "+", 
      #id = td$name, cigar = td$cigar, seq = td$seq, 
      stringsAsFactors = F)
    tr$end = tr$beg + tr$len - 1
    tr$srd[which(td$revstrand)] = "-"
    if(pileup) tr = cbind(tr, y = pileupReads(tr))
    res = rbind(res, tr)
  }
  res
}

rename_genefam <- function(ti) {
  mapping = list(
  	"F-box" = "F-box",
  	"RLK" = c("LRR-RLK", "RLK"),
  	"TE" = "TE",
  	"Unknown" = "Unknown",
  	"NBS-LRR" = c("CC-NBS-LRR", "TIR-NBS-LRR", "NB-ARC", "TIR"),
  	'CRP:NCR' = 'NCR',
  	'CRP:other' = c('CRP0000-1030','CRP1600-6250'),
  	'LRR' = 'LRR',
  	'HSP70' = 'HSP70',
  	'Zinc-Finger' = 'Zinc-Finger',
#  	'Cytochrome' = 'Cytochrome',
  	'Pkinase' = c('Pkinase', 'Pkinase_Tyr')
  )
  to = ti
  for (fam in names(mapping)) {
  	to$fam[to$fam %in% mapping[[fam]]] = fam
  }
  to$fam[! to$fam %in% names(mapping)] = 'Miscellaneous' #'Pfam:other'
  to$fam = factor(to$fam, levels = unique(to$fam))
  to
}

tname = "HM101"
qnames_all = c(
  "HM058", "HM056", "HM125", "HM129", "HM034", 
  "HM095", "HM060", "HM185", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM340", "HM324",
  "HM056.AC", "HM034.AC", "HM340.AC"
)
qnames_12 = c(
  "HM058", "HM056", "HM125", "HM129", "HM034", 
  "HM095", "HM060", "HM185", "HM004", "HM050", 
  "HM023", "HM010"
)
qnames_15 = c(
  "HM058", "HM056", "HM125", "HM129", "HM034", 
  "HM095", "HM060", "HM185", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM340", "HM324"
)
qnames_alpaca = c("HM056.AC", "HM034.AC", "HM340.AC")
qnames_alpaca_comp = c("HM056", "HM056.PJ", "HM056.AC", "HM034", "HM034.PJ", "HM034.AC", "HM340", "HM340.PJ", "HM340.AC", "HM101")
qnames_ingroup = qnames_12
orgs = c(tname, qnames_all)
cfgs = get_genome_cfgs(orgs)
tcfg = cfgs[[tname]]
