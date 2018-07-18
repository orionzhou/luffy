require(GenomicRanges)
require(dplyr)
source('Location.R')

## maize
dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
fg = file.path(dirg, "51.tbl")
tg = read.table(fg, sep = "\t", header = F, as.is = T)
colnames(tg) = c("chr", 'beg', 'end', 'srd', 'id', 'type', 'cat')
tg = tg[tg$chr %in% 1:10,]

grc = reduce(with(tg[tg$type == 'cds',], GRanges(seqnames = chr, ranges = IRanges(beg, end = end))))
gru = reduce(with(tg[tg$type %in% c('utr5','utr3'),], GRanges(seqnames = chr, ranges = IRanges(beg, end = end))))
gru = setdiff(gru, grc)
gri = reduce(with(tg[tg$type == 'intron',], GRanges(seqnames = chr, ranges = IRanges(beg, end = end))))
gri = setdiff(gri, reduce(c(grc, gru)))
stopifnot(sum(width(intersect(grc, gri))) == 0)
stopifnot(sum(width(intersect(grc, gru))) == 0)
stopifnot(sum(width(intersect(gri, gru))) == 0)

flen = file.path(dirg, "15.sizes")
tlen = read.table(flen, sep = "\t", header = F, as.is = T)
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2, stringsAsFactors = F)
tt = tt[tt$chr %in% 1:10,]
grt = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
grx = setdiff(grt, reduce(c(grc, gri, gru)))

twc = data.frame(chr = as.character(seqnames(grc)), beg = start(grc), end = end(grc), type = 'CDS', stringsAsFactors = F)
twu = data.frame(chr = as.character(seqnames(gru)), beg = start(gru), end = end(gru), type = 'UTR', stringsAsFactors = F)
twi = data.frame(chr = as.character(seqnames(gri)), beg = start(gri), end = end(gri), type = 'Intron', stringsAsFactors = F)
twx = data.frame(chr = as.character(seqnames(grx)), beg = start(grx), end = end(grx), type = 'Intergenic', stringsAsFactors = F)
tw = rbind(twc, twi, twu, twx)
tw = tw[order(tw$chr, tw$beg),]
stopifnot(sum(width(grt)) == sum(tw$end-tw$beg+1))

fo = file.path(dirg, "52.itv.tsv")
write.table(tw, fo, sep = "\t", row.names = F, col.names = T, quote = F)

## get longest transcripts
f_gtb = file.path(dirg, "51.gtb")
f_tbl = file.path(dirg, "51.tbl")
tg = read.table(f_gtb, sep = "\t", header = T, as.is = T)[,1:2]
tt = read.table(f_tbl, sep = "\t", header = F, as.is = T)
colnames(tg) = c("tid", "gid")
colnames(tt) = c("chr", "beg", "end", "srd", "tid", "type", "fam")
tt = tt[tt$type %in% c('cds'),]
grp = dplyr::group_by(tt, tid)
tt2 = dplyr::summarise(grp, len = sum(end - beg + 1))

tg2 = merge(tg, tt2, by = 'tid')
grp = dplyr::group_by(tg2, gid)
tg3 = dplyr::summarise(grp, tid = tid[which(len==max(len))[1]])

fo = file.path(dirg, "57.longest.tsv")
write.table(tg3, fo, sep = "\t", row.names = F, col.names = F, quote = F)


## medicago
chrs = sprintf("chr%s", 1:8)

dirr = file.path(Sys.getenv('genome'), "HM101")
fa = file.path(dirr, '15.sizes')
fg = file.path(dirr, '16.gap.bed')
ta = read.table(fa, sep = '\t', header = F, as.is = T)
tg = read.table(fg, sep = '\t', header = F, as.is = T)
ta = ta[ta$V1 %in% chrs,]
tg = tg[tg$V1 %in% chrs,]
ga = GRanges(seqnames = ta$V1, ranges = IRanges(1, end = ta$V2))
gg = GRanges(seqnames = tg$V1, ranges = IRanges(tg$V2+1, end = tg$V3))
gr = setdiff(ga, gg)

fl = file.path(dirr, "51.tbl")
tl = read.table(fl, sep = '\t', header = F, as.is = T)
colnames(tl) = c("chr", "beg", "end", "srd", "id", "type", "cat")
tl = tl[tl$cat != 'TE',]

tc = tl[tl$type == 'cds',]
g_cds = GRanges(seqnames = tc$chr, ranges = IRanges(tc$beg, end = tc$end))
g_cds = reduce(g_cds)
g_cds = intersect(g_cds, gr)

tc = tl[tl$type == 'intron',]
g_ito = GRanges(seqnames = tc$chr, ranges = IRanges(tc$beg, end = tc$end))
g_ito = reduce(g_ito)
g_ito = setdiff(g_ito, g_cds)
g_ito = intersect(g_ito, gr)

tc = tl[tl$type == 'utr5' | tl$type == 'utr3',]
g_utr = GRanges(seqnames = tc$chr, ranges = IRanges(tc$beg, end = tc$end))
g_utr = reduce(g_utr)
g_utr = setdiff(g_utr, union(g_cds, g_ito))
g_utr = intersect(g_utr, gr)

g_ige = setdiff(gr, union(union(g_cds, g_ito), g_utr))


to = data.frame(chr = seqnames(go), beg = start(go) - 1, end = end(go))

types = c('CDS', 'Intron', 'UTR', 'Intergenic')
gos = list(g_cds, g_ito, g_utr, g_ige)

to = data.frame()
for (i in 1:length(types)) {
  type = types[i]
  go = gos[[i]]
  to = rbind(to, 
    data.frame(chr = seqnames(go), beg = start(go), end = end(go), type = type))
}

to = to[order(to$chr, to$beg, to$end), ]
fo = file.path(dirr, "51.merged.tbl")
write.table(to, fo, row.names = F, col.names = T, sep = "\t", quote = F)