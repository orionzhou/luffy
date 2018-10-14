require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(gtable)
#require(Gviz)
source('Location.R')
source('comp.fun.R')

diro = file.path(Sys.getenv("misc3"), "comp.stat")

#qnames = c("HM056", "HM056.AC", "HM340", "HM340.AC")
#qnames = c("HM056", "HM324", "Malbus")

tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

###### basic assembly stats
fr = file.path(diro, "00.sequencing.raw.tbl")
tr = read.table(fr, header = T, sep = "\t", as.is = T)
accs = sapply(strsplit(tr$label, split = "_"), "[", 1)
accs[accs == "HM018"] = "HM023"
idxs = accs %in% qnames
tsa = cbind(org = accs[idxs], tr[idxs,])
tsa$insert_size = as.numeric(tsa$insert_size)

libtype = rep(NA, nrow(tsa))
libtype[! is.na(tsa$insert_size) & tsa$insert_size > 300] = 'LIPE'
libtype[! is.na(tsa$insert_size) & tsa$insert_size <= 300] = 'SIPE'
idxs_sipe = unique(c(grep("SIPE", tsa$label), grep("SIPE", tsa$label.1)))
idxs_lipe = unique(c(grep("LIPE", tsa$label), grep("LIPE", tsa$label.1)))
libtype[is.na(tsa$insert_size) & 1:length(libtype) %in% idxs_sipe] = "SIPE"
libtype[is.na(tsa$insert_size) & 1:length(libtype) %in% idxs_lipe] = "LIPE"

tsa = cbind(tsa, libtype = libtype)
to1 = ddply(tsa, .(org, libtype), summarise, Gb = sum(filtered_count/1000000000*cycles*2))
to = reshape(to1, direction = 'wide', timevar = 'libtype', idvar = 'org')

reflength <- 500000000
to = cbind(to, coverage = (to$Gb.LIPE + to$Gb.SIPE) / (reflength/1000000000))

fr2 = file.path(diro, "00.assembly.raw.tbl")
tr2 = read.table(fr2, header = T, sep = "\t", as.is = T)
tr2$Accession[tr2$Accession == "HM018"] = "HM023"
to = merge(to, tr2[,c('Accession', 'SIPE.status', 'LIPE.status')], by.x = 'org', by.y = 'Accession', all.x = T)
too = data.frame('org' = to$org, 'SIPE' = format(to$Gb.SIPE, digits = 3, nsmall = 1), 'LIPE' = format(to$Gb.LIPE, digits = 3, nsmall = 1), 'SIPE detail' = to$SIPE.status, 'LIPE detail' = to$LIPE.status, 'Total coverage (fold)' = format(to$coverage, digits = 3, nsmall = 1), stringsAsFactors = F)

fo = file.path(diro, "00.stat.tbl")
write.table(too, fo, sep = "\t", row.names = F, col.names = T, quote = F)


library(Biostrings)
source(paste("http://faculty.ucr.edu/~tgirke/",
    "Documents/R_BioCond/My_R_Scripts/contigStats.R", sep=''))

reflength <- 500000000
stats = list()
for (qname in c(tname, qnames)) {
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  flen = file.path(dir, "15.sizes")
  tlen = read.table(flen, sep = "\t", header = F, as.is = T)
  total_len = sum(tlen$V2)
  
  fgap = file.path(dir, "16.gap.bed")
  tgap = read.table(fgap, sep = "\t", header = F, as.is = T)
  total_gap = sum(tgap$V3 - tgap$V2)
  total_bases = total_len - total_gap
  
  if(qname == tname) {
    scf_stat = rep(0, 4)
    ctg_stat = rep(0, 4)
  } else {
    f_scf = file.path(dir, "11_genome.fas")
    assembly <- readDNAStringSet(f_scf, "fasta")
    N <- list(acc = width(assembly))
#    reflength <- sapply(N, sum)
    st <- contigStats(N = N, reflength = reflength, style = "data")
    scf_stat = as.integer(st$Contig_Stats)[c(8,2,6,4)]
    
    f_ctg = file.path(dir, "ctg.fas")
    assembly <- readDNAStringSet(f_ctg, "fasta")
    N <- list(acc = width(assembly))
#    reflength <- sapply(N, sum)
    st <- contigStats(N = N, reflength = reflength, style = "data")
    ctg_stat = as.integer(st$Contig_Stats)[c(8,2,6,4)]
  }
  
  # repeat stats
  frep = file.path(dir, "12.rm.bed")
  trep = read.table(frep, sep = "\t", header = F, as.is = T)
  brep = sum(trep$V3 - trep$V2)
  pct_rep = brep / total_bases * 100
  pct_rep = sprintf("%.02f", pct_rep)

  stat = c(total_len, total_bases, scf_stat, ctg_stat, brep)
  stat = format(stat, big.mark = ",")
  stat = c(stat, pct_rep)
  stats[[qname]] = matrix(stat, nrow = 1, 
    dimnames = list(NULL, c("Total Span", "Total Bases",
    "Number", "N50", "Median", "Max", "Number", "N50", "Median", "Max", 
    "Bases", "Percent")))
}
ds = do.call(rbind.data.frame, stats)
fo = file.path(diro, "01_assembly_stat.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)

for (qname in qnames) {
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  flen = file.path(dir, "15.sizes")
  tlen = read.table(flen, sep = "\t", header = F, as.is = T)
  total_len = sum(tlen$V2)
  cat(qname, sum(tlen$V2[tlen$V2 >= 100000]) / total_len, "\n")
}

tc = read.table(tcfg$gene, header = F, sep = "\t", as.is = T)
colnames(tc) = c("chr", "beg", "end", "srd", "id", "type", "fam")
grc = with(tc[tc$type == 'cds',], GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
stats = list()
for (qname in qnames) {
  fgax = sprintf("%s/%s_HM101/23_blat/31.5/gax", Sys.getenv("misc3"), qname)
  tgax = read.table(fgax, header = F, sep = "\t", as.is = T)[,1:3]
  grg = with(tgax, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
  
  bp_cov = sum(width(reduce(grg)))
  bp_cov_cds = sum(width(GenomicRanges::intersect(grg, grc)))
  
  stats[[qname]] = matrix(c(bp_cov, bp_cov_cds), nrow = 1, dimnames = list(NULL,
    c("bp_cov", "bp_cov_cds")))
}
ds = do.call(rbind.data.frame, stats)
ungap = sum(tt$end) - sum(tp$end - tp$beg + 1)
ds = cbind(ds, pct_cov = ds$bp_cov / ungap, pct_cov_cds = ds$bp_cov_cds / sum(width(reduce(grc))))

##### functional annotation stat
stats = list()
statd = list()
nstats = list()
rstats = list()
cstats = list()
tstats = list()

orgs_rnaseq = c("HM034", "HM056", "HM101", "HM340")

fi = file.path(Sys.getenv("misc2"), "genefam", "nbs.info")
nfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("misc2"), "genefam", "crp.info")
cfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("misc2"), "genefam", "rlk.info")
rfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("data"), 'db', 'pfam', 'genefam.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
tfams = ti$dom[ti$fam == 'TE']
fams = unique(ti$fam[order(ti$pri)])

for (qname in c(tname, qnames_15)) {
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  fg = file.path(dir, "51.gtb")
  tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,15:17)]
  ngene = nrow(tg)
  ids_nonte = tg$id[tg$cat2 != "TE"]
  n_nonte = length(ids_nonte)
  
  dtn = table(tg$cat2)
  x = c()
  x[fams] = 0
  x[names(dtn)] = dtn
  statd[[qname]] = matrix(c(x, ngene),
    nrow = 1, dimnames = list(NULL, c(fams, "Total Genes")))
  
  tf = file.path(dir, "51.fas")
  y = read.fasta(tf, as.string = TRUE, seqtype = "AA")
  dl = ldply(y, nchar)
  lens = dl$V1[dl$.id %in% ids_nonte]
  mean_prot = mean(lens)
  med_prot = median(lens)

  tt = tg[tg$cat2 == 'TE',]
  dtn = table(tt$cat3)
  dtn = dtn[names(dtn) != '']
  x = c()
  x[tfams] = 0
  x[names(dtn)] = dtn
  tstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, tfams))
  n_te = nrow(tt)
  
  tn = tg[tg$cat2 %in% c('CC-NBS-LRR', 'TIR-NBS-LRR'),]
  dtn = table(tn$cat3)
  x = c()
  x[nfams] = 0
  x[names(dtn)] = dtn
  nstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, nfams))
  n_nbs = nrow(tn)
  
  tr = tg[tg$cat2 %in% c('RLK','LRR-RLK'),]
  dtn = table(tr$cat3)
  x = c()
  x[rfams] = 0
  x[names(dtn)] = dtn
  rstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, rfams))
  n_rlk = nrow(tr)
  
  tc = tg[tg$cat2 %in% c('CRP0000-1030','NCR','CRP1600-6250'),]
  dtn = table(tc$cat3)
  x = c()
  x[cfams] = 0
  x[names(dtn)] = dtn
  cstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, cfams))
  n_ncr = sum(tc$cat2 == 'NCR')
  
  tf = tg[tg$cat2 %in% c('F-box'),]
  n_fbx = nrow(tf)
  
  ids = ids_nonte
  ids_rna = c(); ids_hom = c()
  if(qname %in% orgs_rnaseq) {
    org = qname
    ff = file.path(Sys.getenv("misc2"), "rnaseq/mt/31_cufflinks", org, 'isoforms.fpkm_tracking')
    tf = read.table(ff, header = T, sep = "\t", as.is = T)
    ids_rna = tf$tracking_id[tf$tracking_id %in% ids & tf$FPKM > 0]
    p_rna = sprintf("%.01f", length(ids_rna) / length(ids) * 100)
  } else {
    p_rna = '-'
  }
  if(qname != tname) {
    ft = sprintf("%s/%s_HM101/51_ortho/31.ortho.tbl", Sys.getenv("misc3"), qname)
    tt = read.table(ft, header = T, sep = "\t", as.is = T)
    ids_hom = tt$qid[tt$qid %in% ids]
    p_hom = sprintf("%.01f", length(ids_hom) / length(ids) * 100)
  } else {
    p_hom = '-'
  }
  p_sup = sprintf("%.01f", length(unique(c(ids_rna, ids_hom))) / length(ids) * 100)
  
  stats[[qname]] = matrix(c(ngene, n_te, n_nonte, n_nbs, n_fbx, n_rlk, n_ncr, med_prot, p_rna, p_hom, p_sup), nrow = 1, dimnames = list(NULL, c("# Total Genes", 'TE', 'non-TE','NBS-LRR','F-box','LRR-RLK','NCR',"Median Prot Length", 'RNA-seq (%)', 'Homology (%)', 'RNA-seq + Homology (%)')))
}
ds = do.call(rbind.data.frame, stats)
#for (i in 1:ncol(do)) { do[,i] = format(do[,i], big.mark = ",") }
fo = file.path(diro, "03_annotation.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)

ds = do.call(rbind.data.frame, statd)
fo = file.path(diro, "04.annotation.detail.tbl")
write.table(t(ds), fo, sep = "\t", row.names = T, col.names = T, quote = F)

ds = do.call(rbind.data.frame, nstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "04.nbs.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

ds = do.call(rbind.data.frame, rstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "04.rlk.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

ds = do.call(rbind.data.frame, cstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "04.crp.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

ds = do.call(rbind.data.frame, tstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "04.te.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

##### comparative / novel-seq stats
stats = list()
for (qname in qnames) {
  #add repeatmasker stats
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  flen = file.path(dir, "15.sizes")
  tlen = read.table(flen, sep = "\t", header = F, as.is = T)
  total_len = sum(tlen$V2)
  
  fgap = file.path(dir, "16.gap.bed")
  tgap = read.table(fgap, sep = "\t", header = F, as.is = T)
  total_gap = sum(tgap$V3 - tgap$V2)
  total_bases = total_len - total_gap
  
  fcds = file.path(dir, "51.tbl")
  tcds = read.table(fcds, sep = "\t", header = F, as.is = T)[,1:6]
  tcds = tcds[tcds$V6 == 'cds',]
  grc = GRanges(seqnames = tcds$V1, ranges = IRanges(tcds$V2, end = tcds$V3))
  grc = reduce(grc)
  
  frep = file.path(dir, "12.rm.bed")
  trep = read.table(frep, sep = "\t", header = F, as.is = T)
  brep = sum(trep$V3 - trep$V2)
  pct_rep = brep / total_bases * 100

  dir = sprintf("%s/%s_%s/23_blat", Sys.getenv("misc3"), qname, tname)
  fi = file.path(dir, '41.5/gax')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:3,10)]
  colnames(ti) = c('tid', 'tbeg', 'tend', 'lev')
  aligned = sum(ti$tend - ti$tbeg + 1)
#  pct_aligned = aligned / total_bases * 100

  fi = file.path(dir, '31.9/gax')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:3,10)]
  colnames(ti) = c('tid', 'tbeg', 'tend', 'lev')
  ti = ti[ti$lev <= 2,]
  synteny = sum(ti$tend - ti$tbeg + 1)
#  pct_synteny = synteny / total_bases * 100

  dir = sprintf("%s/%s_%s/41_novseq", Sys.getenv("misc3"), qname, tname)
  fi = file.path(dir, '21.bed')
  ti = read.table(fi, sep = "\t", header = F, as.is = T)
  bnov = sum(ti$V3 - ti$V2)
  pnov = bnov / total_bases * 100
  
  grn = GRanges(seqnames = ti$V1, ranges = IRanges(ti$V2+1, end = ti$V3))
  grn = reduce(grn)
  bcds = sum(width(intersect(grc, grn)))
  pcds = bcds / bnov * 100

  total_bases = format(total_bases, big.mark = ",")
  brep = format(brep, big.mark = ",")
#  pct_rep = sprintf("%.01f%%", pct_rep)
  aligned = format(aligned, big.mark = ",")
#  pct_aligned = sprintf("%.01f%%", pct_aligned)
  synteny = format(synteny, big.mark = ",")
#  pct_synteny = sprintf("%.01f%%", pct_synteny)
  bnov = format(bnov, big.mark = ",")
  pnov = sprintf("%.01f%%", pnov)
  bcds = format(bcds, big.mark = ",")
  pcds = sprintf("%.01f%%", pcds)
  stats[[qname]] = matrix(c(total_bases, brep, aligned, synteny, bnov, pnov, bcds, pcds), 
    nrow = 1, dimnames = list(NULL, c("Total Bases", "Repetitive", "Alignable to HM101", "Bases in Synteny Blocks", "Novel Sequences", "", "Novel Coding Seq", "")))
}
ds = do.call(rbind.data.frame, stats)
fo = file.path(diro, "11_comp_stat.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)

##### variant stats
stats_bp = list()
stats_cnt = list()
stats_snpd = list()
for (qname in qnames) {
  dir = sprintf("%s/%s_%s", Sys.getenv("misc3"), qname, tname)
  fi = file.path(dir, '23_blat/31.9/gax')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:3,10)]
  colnames(ti) = c('tid', 'tbeg', 'tend', 'lev')
  ti = ti[ti$lev == 1,]
  ali = sum(ti$tend - ti$tbeg + 1)
  
  fi = file.path(dir, '23_blat/31.9/snp')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:2,8)]
  colnames(ti) = c('chr', 'beg', 'lev')
  snp = sum(ti$lev == 1)
  snpd = snp / ali

  fv = file.path(dir, '31_sv/05.stx')
  tv = read.table(fv, header = T, sep = "\t", as.is = T)

  ti = tv[tv$type == 'INS',]
  ti = cbind(ti, len = ti$qend - ti$qbeg + 1)
  tis = ti[ti$len < 50,]
  si_n = nrow(tis); si_b = sum(tis$len)
  til = ti[ti$len >= 50,]
  li_n = nrow(til); li_b = sum(til$len)

  ti = tv[tv$type == 'CNG',]
  ti = cbind(ti, len = ti$qend - ti$qbeg + 1)
  g_n = nrow(ti); g_b = sum(ti$len)

  ti = tv[tv$type == 'DEL',]
  ti = cbind(ti, len = ti$tend - ti$tbeg + 1)
  tis = ti[ti$len < 50,]
  sd_n = nrow(tis); sd_b = sum(tis$len)
  til = ti[ti$len >= 50,]
  ld_n = nrow(til); ld_b = sum(til$len)

  ti = tv[tv$type == 'CNL',]
  ti = cbind(ti, len = ti$tend - ti$tbeg + 1)
  l_n = nrow(ti); l_b = sum(ti$len)
  
  ti = tv[tv$type %in% c('TLC:INS', 'TLC:DEL'),]
  gri = with(ti, GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend)))
  grr = reduce(gri)
  t_n = length(grr); t_b = sum(width(grr))
  
#  stat = c(snp, si_n, si_b, sd_n, sd_b, li_n, li_b, ld_n, ld_b,
#    g_n, g_b, l_n, l_b, t_n, t_b)
#  stat = format(stat, big.mark = ",")
#  stat = c(stat[1], snpd, stat[-1])
  cnames = c("SNP", "sIns", "sDel", "lIns", "lDel", "CNG", "CNL", "TLC")
  stats_cnt[[qname]] = matrix(c(snp, si_n, sd_n, li_n, ld_n, g_n, l_n, t_n), nrow = 1, dimnames = list(NULL, cnames))
  stats_bp[[qname]] = matrix(c(snp, si_b, sd_b, li_b, ld_b, g_b, l_b, t_b), nrow = 1, dimnames = list(NULL, cnames))
  stats_snpd[[qname]] = snpd
}
tcn = do.call(rbind.data.frame, stats_cnt)
tbp = do.call(rbind.data.frame, stats_bp)
tsd = data.frame(org = qnames, snpd = as.numeric(unlist(stats_snpd[qnames])))

t1 = cbind(org = sprintf("%s(#)", qnames), format(tcn, big.mark = ","))
t2 = cbind(org = sprintf("%s(bp)", qnames), format(tbp, big.mark = ","))
t2$SNP = ''
t3 = rbind(t1, t2)

snpd = rep('', nrow(t3))
snpd[t3$SNP != ''] = sub("^(-?)0.", "\\1.", sprintf("%.04f", tsd$snpd))
t4 = cbind(t3[,1:2], snpd = snpd, t3[,3:ncol(t3)])
t4$snpd[t4$SNP==''] = ''
colnames(t4)[c(3:7,10)] = c("SNP Density", "Small Ins", "Small Del", "Large Ins", "Large Del", "Translocation")

fo = file.path(diro, "15_vnt_stat.tbl")
write.table(t4, fo, sep = "\t", row.names = F, col.names = T, quote = F)

##### correlation of genome size estimates with Kelly's data
tsa = data.frame(org = c(tname, qnames), tspan = NA, tbase = NA, stringsAsFactors = F)
for (i in 1:nrow(tsa)) {
  org = tsa$org[i]
  dir = sprintf("%s/%s", Sys.getenv("genome"), org)
  flen = file.path(dir, "15.sizes")
  tlen = read.table(flen, sep = "\t", header = F, as.is = T)
  tsa$tspan[i] = sum(tlen$V2)
  
  fgap = file.path(dir, "16.gap.bed")
  tgap = read.table(fgap, sep = "\t", header = F, as.is = T)
  total_gap = sum(tgap$V3 - tgap$V2)
  tsa$tbase[i] = tsa$tspan[i] - total_gap
}

fsk = file.path(diro, "genome_size", "01.kelly.tbl")
tsk1 = read.table(fsk, header = T, sep = "\t", as.is = T)
tsk = ddply(tsk1, .(org), summarise, wt = mean(wt))
tsk$org[tsk$org == 'HM029'] = 'HM340'
tsk$org[tsk$org == 'HM018'] = 'HM023'

to = merge(tsa, tsk, by = "org")
to = cbind(to, tspanr = to$tspan / mean(to$tspan))
fit = lm(tspanr ~ wt, data = to)
fitsum = summary(fit)
f = fitsum$fstatistic
pval = pf(f[1], f[2], f[3], lower.tail = F)
attributes(pval) <- NULL
txt = sprintf("'R'^2*' = %.03f  p-value = %.04f'", fitsum$r.squared, pval)

corc = cor(to$tspanr, to$wt)
pval = cor.test(to$tspanr, to$wt)$p.value
txt = sprintf("correlation coefficient = %.03f  p-value = %.04f", corc, pval)

p = ggplot(to) +
  geom_point(aes(x = wt, y = tspanr), shape = 4, size = 2) +
#  geom_text(aes(x = n_org, y = 0, label = org), geom_params=list(size = 2.5, vjust = 0, angle = 30)) +
  stat_smooth(aes(x = wt, y = tspanr), method = 'lm', fill = 'azure3', size = 0.5) +
#  scale_shape(name = "", solid = FALSE, guide = F) +
#  scale_color_manual(name = "", labels = c('Core-genome', 'Pan-genome'), values = c("dodgerblue", "firebrick1"), guide = guide_legend(label.position = "left", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  scale_x_continuous(name = 'Genome size estimate (fluorometry)') +
  scale_y_continuous(name = 'AllPaths assembly size (normailized)', expand = c(0.02, 0)) + 
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.2, 'lines')) +
  theme(legend.position = "top", legend.key.size = unit(1, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,1,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 10, angle = 0)) +
  theme(axis.title.y = element_text(size = 10, angle = 90)) +
  theme(axis.text.x = element_text(size = 9, color = "dodgerblue4")) +
  theme(axis.text.y = element_text(size = 9, color = "dodgerblue4", angle = 0, hjust  = 0.5)) +
  annotate("text", x = 0.97, y = 1.05, size = 3, hjust = 0, label = txt, parse = F)

fp = file.path(diro, "genome.size.pdf")
ggsave(p, filename = fp, width = 6, height = 5)

