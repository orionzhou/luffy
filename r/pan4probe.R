require(GenomicRanges)
require(plyr)
require(rtracklayer)

dir = '/home/youngn/zhoup/Data/misc3/pan4probe'

dirs = '/home/youngn/zhoup/Data/genome/pan4/18_stat_k60'
fbwg = file.path(dirs, "11_gc.bw")
fbwm = file.path(dirs, "15_mapp.bw")
fbwd = file.path(dirs, "17_deltag.bw")

### partition genomic regions by gc/mp/dg
com = matrix( c(
  0.3, 0.583, 1.000, 1.000, -7.695, 6.47,
  0.3, 0.583, 0.500, 0.500, -7.695, 6.47,
  0.3, 0.583, 0.333, 0.333, -7.695, 6.47,
  0.3, 0.583, 0.250, 0.250, -7.695, 6.47,
  0.3, 0.583, 0.200, 0.200, -7.695, 6.47,
  0.3, 0.583, 0.001, 0.199, -7.695, 6.47
  ), ncol = 6, byrow = T)
co = data.frame(com)
colnames(co) = c('gc_min', 'gc_max', 'mp_min', 'mp_max', 'dg_min', 'dg_max')

fp = file.path(dirs, "../17.5.gap.bed")
tp = read.table(fp, header = F, sep = "\t", as.is = T, quote = "")
colnames(tp) = c("chr", 'beg', 'end')
gp = GRanges(seqnames = Rle(tp$chr), ranges = IRanges(tp$beg + 1, end = tp$end))

ff = file.path(dir, "15.merged.bed")
tf = read.table(ff, header = F, sep = "\t", as.is = T)
colnames(tf) = c("chr", "beg", "end")
gf = GRanges(seqnames = tf$chr, ranges = IRanges(tf$beg + 1, end = tf$end))

gr = gf
bwg = import(fbwg, which = gr, asRangedData = F)
bwm = import(fbwm, which = gr, asRangedData = F)
bwd = import(fbwd, which = gr, asRangedData = F)
chrs = seqnames(bwg)
poss = start(bwg)
gc = bwg$score
rm(bwg)
mp = bwm$score
rm(bwm)
dg = bwd$score
rm(bwd)

err = 1e-4
for (i in 1:nrow(co)) {
  idxs = 
    gc > co$gc_min[i] - err & gc < co$gc_max[i] + err &
    mp > co$mp_min[i] - err & mp < co$mp_max[i] + err &
    dg > co$dg_min[i] - err & dg < co$dg_max[i] + err

  gr = GRanges(seqnames = chrs[idxs], ranges = 
    IRanges(poss[idxs], end = poss[idxs]))
  grr = setdiff(reduce(gr), gp)

  tr = data.frame(chr = as.character(seqnames(grr)), 
    beg = start(grr) - 1, end = end(grr))
  fo = sprintf("%s/21.cds.nov/rd%d.bed", dir, i)
  options(scipen = 999)
  write.table(tr, fo, col.names = F, row.names = F, sep='\t', quote = F)
}

### probe stat per gene
probe_stat <- function(fi) {
  t = read.table(fi, sep = "\t", header = T, as.is = T)[, c(1:5,7)]
  t$prbnum[t$prbnum > 3] = 3
  table(t[, c('type', 'prbnum')])
}
probe_stat(sprintf("%s/21.cds.nov/rd%d/31.tbl", dir, 4))
probe_stat(sprintf("%s/41.tbl", dir))

fp = sprintf("%s/21.cds.nov/rd%d/31.tbl", dir, 4)
t = read.table(fp, sep = "\t", header = T, as.is = T)[, c(1:5,7)]
ts = t[t$type == 'CRP' & t$prbnum == 0,]
sprintf("%s:%d-%d", ts$chr, ts$beg, ts$end)

### collect probe statistics
fp = file.path(dir, "36.prb")
fo = file.path(dir, "37.prbstat.tbl")
fp = file.path(dir, "06.prb")
fo = file.path(dir, "07.prbstat.tbl")

tp = read.table(fp, header = F, sep="\t", as.is=T, quote="")
colnames(tp) = c("chr", "beg")
gr = GRanges(seqnames = Rle(tp$chr), ranges = IRanges(tp$beg, end = tp$beg))

bwg = import(fbwg, which = gr, asRangedData = F)
bwm = import(fbwm, which = gr, asRangedData = F)
bwd = import(fbwd, which = gr, asRangedData = F)
gc = format(bwg$score, digits = 3)
rm(bwg)
mp = format(bwm$score, digits = 3)
rm(bwm)
dg = format(bwd$score, digits = 3)
rm(bwd)

id = sprintf("mtpan4_%06d", 1:nrow(tp))
to = cbind(tp, end = tp$beg + 59, id = id, len = 60, gc = gc, mp = mp, dg = dg)
options(scipen = 999)
write.table(to, fo, col.names = T, row.names = F, sep='\t', quote = F)


### design pan4 back-bone probes
# seqgap.pl -i 11_genome.fas -o 17.1.gap.bed -m 1
# flankBed -l 59 -r 0 -i 17.1.gap.bed > 17.2.gapflank.bed
# cat 17.1.gap.bed 17.2.gapflank.bed | sortBed -i stdin | mergeBed -i stdin > 17.5.gap.bed
# subtractBed -a 15.bed -b 17.5.gap.bed > 17.7.nogap.bed

fa = file.path(dirs, "../17.7.nogap.bed")
ta = read.table(fa, header = F, sep = "\t", as.is = T, quote = "")
colnames(ta) = c("chr", 'beg', 'end')
ga = GRanges(seqnames = Rle(ta$chr), ranges = IRanges(ta$beg + 1, end = ta$end))

gr = ga

bwg = import(fbwg, which = gr, asRangedData = F)
chrs = seqnames(bwg)
poss = start(bwg)
gc = bwg$score
gc_min = 0.3
gc_max = 0.583
idxs_gc = gc >= gc_min & gc <= gc_max
rm(bwg, gc)


bwm = import(fbwm, which = gr, asRangedData = F)
mapp = bwm$score
mapp_min = 1
mapp_max = 1
idxs_mapp = mapp >= mapp_min & mapp <= mapp_max
rm(bwm, mapp)


bwd = import(fbwd, which = gr, asRangedData = F)
deltag = bwd$score
deltagf = deltag[deltag < 900]
length(deltagf)
summary(deltagf)
deltag_min = mean(deltagf) - 2*sd(deltagf)
deltag_max = max(deltagf)
idxs_deltag = deltag >= deltag_min & deltag <= deltag_max
rm(bwd, deltag, deltagf)

idxs = idxs_gc & idxs_mapp & idxs_deltag

gs = GRanges(seqnames = chrs[idxs], ranges = IRanges(poss[idxs], 
  end = poss[idxs]))
gsr = reduce(gs)

tps = data.frame(chr = as.character(seqnames(gsr)), beg = start(gsr) - 1, 
  end = end(gsr))

fuq = file.path(dir, "03.uniq.bed")
options(scipen = 999)
write.table(tps, fuq, col.names = F, row.names = F, sep='\t', quote = F)

