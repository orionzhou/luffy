require(plyr)
require(dplyr)
require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(RColorBrewer)
source("Location.R")
source("comp.fun.R")

dirw = file.path(Sys.getenv("misc3"), "comp.sv.valid")

##### filter SVs (>= 50bp)
org = "HM034"
for (org in c("HM034", "HM056", "HM340")) {
dirg = file.path(Sys.getenv("genome"), org)

dirv = sprintf("%s/%s_HM101/31_sv", Sys.getenv("misc3"), org)
fv = file.path(dirv, "05.stb")
tv = read.table(fv, header = T, sep = "\t", as.is = T)
tv = tv[tv$tlen >= 50 | tv$qlen >= 50,]

fo = sprintf("%s/01_sv/%s.tbl", dirw, org)
write.table(tv, fo, sep = "\t", row.names = F, col.names = T, quote = F)
}

##### classify SVs by gene/TE/intergenic, or SV types
org = "HM034"
for (org in c("HM034", "HM056", "HM340")) {
	fi = sprintf("%s/01_sv/%s.tbl", dirw, org)
	ti = read.table(fi, header = T, sep = "\t", as.is = T)
	
	types = c()
	for (i in 1:nrow(ti)) {
		tlen = ti$tlen[i]; tinfo = ti$tinfo[i]; qlen = ti$qlen[i]; qinfo = ti$qinfo[i]
		is_t = tlen > 0 & tlen > qlen
		info = ifelse(is_t, tinfo, qinfo)
		arys = strsplit(info, split = ",")[[1]]
		x = strsplit(arys, split = "-")
		stypes = sapply(x, "[", 3)
		slens = as.numeric(sapply(x, "[", 2)) - as.numeric(sapply(x, "[", 1)) + 1
		plen = sum(slens[stypes == 'pav'])
		clen = sum(slens[stypes == 'cnv'])
		is_p = plen > clen
		tag = ifelse(is_t, ifelse(is_p, 'del', 'cnl'), ifelse(is_p, 'ins', 'cng'))
		types = c(types, tag)
	}
	table(types)
	
	dirg = file.path('/home/youngn/zhoux379/data/genome', 'HM101')
	fg = sprintf("%s/51.tbl", dirg)
	tg = read.table(fg, header = F, sep = "\t", as.is = T)
	tg = tg[tg$V6 == 'mrna',]
	grg = with(tg, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
	grsv = with(ti, GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend)))
	idxs = intersect_idx_maxovlp(grsv, grg)
	
	gtags = rep('Intergenic', nrow(ti))
	idxs_sv = which(!is.na(idxs))
	gtags[idxs_sv] = tg$V7[idxs_sv]
	gtags[!gtags %in% c("Intergenic", "TE", "Unknown")] = 'Gene'
	table(gtags)
	
	to = cbind(ti, stype = types, gtype = gtags)
	fo = sprintf("%s/01_sv/%s.1.tsv", dirw, org)
	write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
}

##### read in pacbio BAM for support (obsolete - see comp.sv.valid.py)
f_bam = sprintf("%s/pacbio/%s_%s/15.bam", Sys.getenv("misc3"), org, org)
bam = bamReader(f_bam, idx = T)

to = cbind(tv, nr1 = NA, nr2 = NA)
for (i in 1:100) {
#  i = 1
  chr = tv$qchr[i]; beg = tv$qbeg[i]; end = tv$qend[i]
  b1 = max(1, beg - 10)
  e1 = beg + 10
  b2 = max(1, end - 10)
  e2 = end + 10
  gr1 = GRanges(seqnames = chr, ranges = IRanges(b1, end = e1))
  gr2 = GRanges(seqnames = chr, ranges = IRanges(b2, end = e2))
  
  ds1 = read_bam(bam, gr1)
  ds2 = read_bam(bam, gr2)
  to$nr1[i] = nrow(ds1)
  to$nr2[i] = nrow(ds2)
}

sum(to$nr1 > 1 & to$nr2 > 1)
sum(to$nr1 <= 1 & to$nr2 <= 1)

#### generate stat table
txt = paste("", "Deletion", "Insertion", "Copy Number Loss", "Copy Number Gain", "Gene", "TE", "Unknown", "Intergeic", sep = "\t")
txts = c(txt)
for (org in c("HM034", "HM056", "HM340")) {
  fi = sprintf("%s/01_sv/%s.1.tsv", dirw, org)
  ti = read.table(fi, header = T, sep = "\t", as.is = T)
  fp = sprintf("%s/02_pacbio/%s.tbl", dirw, org)
  tp = read.table(fp, header = T, sep = "\t", as.is = T)
  vtag = tp$n1 >= 5 & tp$n2 >= 5
  tt = cbind(ti, vtag = vtag)
  ds1 = ddply(tt, .(stype), summarise, n = length(vtag), nv = sum(vtag)/n)
  ds2 = ddply(tt, .(gtype), summarise, n = length(vtag), nv = sum(vtag)/n)
  ds1 = ds1[match(ds1$stype, c("del", "ins", "cnl", "cng")),]
  ds2 = ds2[match(ds2$gtype, c("Gene", "TE", "Unknown", "Intergenic")),]
  txt = paste(c(org, 
  	sprintf("%.01f%% (%s)", ds1$nv*100, prettyNum(ds1$n, big.mark = ",")),
  	sprintf("%.01f%% (%s)", ds2$nv*100, prettyNum(ds2$n, big.mark = ","))),
  sep = "\t", collapse = "\t")
  txts = c(txts, txt)
}
write(txts, file = file.path(dirw, '11.stat.tsv'))
