require(dplyr)
require(tidyr)
require(GenomicRanges)

dirw = '/home/springer/zhoux379/data/genome/B73'

## gene longest boundary
fi = file.path(dirw, "v37/t2.gtb")
ti = read.table(fi, sep = "\t", as.is = T, header = T)
colnames(ti)[2] = 'gid'

grp = dplyr::group_by(ti, gid)
tg = dplyr::summarise(grp, chr = chr[1], beg = min(beg), end = max(end), srd = srd[1])

to = tg[,c(2:4,1)]
to$beg = to$beg - 1
to = to[order(to$chr, to$beg),]

fo = file.path(dirw, "v37/gene.bed")
write.table(to, fo, sep = "\t", row.names = F, col.names = F, quote = F)

## gene length for HT-Seq
fi = file.path('/home/springer/zhoux379/data/genome/B73', "51.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)
ti = ti[ti$etype == 'exon',]

#gr = with(ti, GRanges(seqnames = seqid, ranges = IRanges(beg, end = end), gid = gid))
#x = unlist(reduce(split(gr, elementMetadata(gr)$gid)))
#tg = data.frame(gid = names(x), chr = seqnames(x), beg = start(x), end = end(x),  stringsAsFactors = F) %>% 
tg = ti %>% 
	group_by(gid, tid) %>% 
	summarise(size = sum(end - beg + 1)) %>%
	group_by(gid) %>%
	summarise(size = max(size))

## gap stats
orgs = c("B73", "W22", "PH207")

tp = data.frame()
for (org in orgs) {
	fi = sprintf("/home/springer/zhoux379/data/genome/%s/16.gap.bed", org)
	ti = read.table(fi, header = F, sep = "\t", as.is = T)
	colnames(ti) = c("chr", "beg", "end")
	tp = rbind(tp, cbind(org = org, ti))
}
tp = cbind(tp, size = tp$end - tp$beg)

grp = group_by(tp, org)
tps = summarise(grp, num = n(), mean = mean(size), median = median(size), min = min(size), max = max(size), total = sum(size))