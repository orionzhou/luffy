require(plyr)
require(dplyr)
require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(RColorBrewer)
source("comp.fun.R")

dirw = file.path(Sys.getenv('misc3'), 'comp.panseq')

##### refine mugsy result
fi = file.path(dirw, '22.tbl')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
#ddply(ti, .(org), summarise, len = sum(end - beg + 1))

to = data.frame()
cid_max = max(ti$cid)
for (qname in qnames) {
#qname = "HM004"
idxs1 = which(ti$org == qname)
t1 = ti[idxs1,]
#t1 = t1[order(t1$chr, t1$beg),]

gr1 = with(t1, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
#grl1 = with(t1, makeGRangesListFromFeatureFragments(seqnames = chr, 
#  fragmentStarts = lapply(beg, rep), fragmentEnds = lapply(end, rep), 
#  strand = rep("*", nrow(t1))))

c1 = coverage(gr1)
l1 = mapply(function(x) {
  matrix(c(start(x), end(x), runValue(x)), nrow = 3, byrow = T)
}, c1)
lens = mapply(ncol, l1)
chrids = rep(names(l1), lens)
m1 = matrix(unlist(l1, use.names = F), ncol = 3, byrow = T)
colnames(m1) = c('beg', 'end', 'val')
t2 = data.frame(m1, stringsAsFactors = F)
t2 = cbind(chr = chrids, t2)
t3 = t2[t2$val > 1,]

idxs_rm = c()
for (i in 1:nrow(t3)) {
  chr = t3$chr[i]; beg = t3$beg[i]; end = t3$end[i]
  idxs = which(t1$org == qname & t1$chr == chr & t1$beg <= end & t1$end >= beg)
  stopifnot(length(idxs) > 1)
  lens = t1$end[idxs] - t1$beg[idxs] + 1
  idx_max = which(lens == max(lens))[1]
  idxs_rm = c(idxs_rm, idxs[-idx_max])
}
t2 = t1[-idxs_rm, ]
gr2 = with(t2, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
stopifnot(sum(width(gr2)) == sum(width(reduce(gr2))))

fb = sprintf("%s/%s_%s/41_novseq/21.bed", Sys.getenv("misc3"), qname, tname)
tb = read.table(fb, sep = "\t", header = F, as.is = T)
colnames(tb) = c('chr', 'beg', 'end', 'type')
tb = within(tb, {beg = beg + 1; chr = sprintf("%s_%d_%d", chr, beg, end);
  end = end - beg + 1; beg = 1})
grb = with(tb, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

grr = setdiff(grb, gr2)
t2c = data.frame(cid = cid_max+1:length(grr), alen = width(grr), org = qname,
  chr = seqnames(grr), beg = start(grr), end = end(grr), srd = "+", 
  stringsAsFactors = F)
t4 = rbind(t2, t2c)
gr4 = with(t4, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
stopifnot(sum(width(gr4)) == sum(width(grb)))
cat(qname, ": ", sum(width(grr)), " overlap removed\n", sep = '')
to = rbind(to, t4)
}

fo = file.path(dirw, "31.refined.tbl")
write.table(to, file = fo, sep = "\t", row.names = F, col.names = T, quote = F)

