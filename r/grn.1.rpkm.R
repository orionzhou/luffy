require(dplyr)
require(GenomicRanges)

dirw = file.path(Sys.getenv("misc2"), "grn23")

### compute RPM
fi = file.path(dirw, '00.5.htseq.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

gids = c()
for (i in 1:nrow(ti)) {
	fh = ti$HtseqFile[i]
	th = read.table(fh, header = F, sep = "\t", as.is = T)
	ngene = nrow(th) - 5
	gids0 = th$V1[-c(ngene+1:5)]
	vals = th$V2[-c(ngene+1:5)]
	vals = vals / (sum(vals) / 1000000)
	if(i == 1) {
		gids = gids0
		to = data.frame(gid = gids, vals = vals, stringsAsFactors = F)
	} else {
		stopifnot(identical(gids, gids0))
		to = cbind(to, vals = vals)
	}
	colnames(to)[i+1] = ti$SampleID[i]
}

dim(to)
fo = file.path(dirw, '33.rpm.tsv')
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### compute RPKM
dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
f_gtb = file.path(dirg, "51.gtb")
f_tbl = file.path(dirg, "51.tbl")
tg = read.table(f_gtb, sep = "\t", header = T, as.is = T)[,1:2]
tt = read.table(f_tbl, sep = "\t", header = F, as.is = T)
colnames(tg) = c("tid", "gid")
colnames(tt) = c("chr", "beg", "end", "srd", "tid", "type", "fam")
tt2 = tt[tt$type %in% c('cds', 'utr5', 'utr3'),]
tg2 = merge(tg, tt2, by = 'tid')

gr = with(tg2, GRanges(seqnames = chr, ranges = IRanges(beg, end = end), gid = gid))
x = unlist(reduce(split(gr, elementMetadata(gr)$gid)))
tr = data.frame(gid = names(x), chr = seqnames(x), beg = start(x), end = end(x), stringsAsFactors = F)
grp = dplyr::group_by(tr, gid)
tr2 = dplyr::summarise(grp, len = sum(end - beg + 1))

to2 = merge(to, tr2, by = 'gid')
stopifnot(nrow(to) == nrow(to2))
for (i in 2:(ncol(to2)-1)) {
	to2[,i] = to2[,i] / (to2$len/1000)
}
to3 = to2[,-ncol(to2)]
dim(to3)
fo = file.path(dirw, '34.rpkm.tsv')
write.table(to3, fo, sep = "\t", row.names = F, col.names = T, quote = F)

