require(plyr)
require(dplyr)

dirw = '/home/springer/zhoux379/data/misc1/rosa/applegenes'

### get contig positions on chrs
fi = file.path(dirw, '31.ctg.gff')
ti = read.table(fi, sep = "\t", as.is = T, header = F)
x = regexec("ID=([\\w\\.]+);", ti$V9, perl = T)
begs = sapply(x, "[", 2)
lens = sapply(x, myfunc <- function(x) attr(x, "match.length")[2])
ids = substr(ti$V9, begs, begs + lens - 1)

x = regexec("Parent=([\\w\\.]+)", ti$V9, perl = T)
begs = sapply(x, "[", 2)
lens = sapply(x, myfunc <- function(x) attr(x, "match.length")[2])
pas = substr(ti$V9, begs, begs + lens - 1)

ti2 = cbind.data.frame(ti[,c(1,3:5,7)], id = ids, pa = pas, stringsAsFactors = F)
colnames(ti2) = c("chr", "type", "beg", "end", "srd", "cid", "pid")

### get primary-ht scaffolds
fp = file.path(dirw, 'Malus_x_domestica.v3.0.a1_pht.tsv')
tp = read.table(fp, sep = "\t", as.is = T, header = F)

###
fa = file.path(dirw, "gene.aln.gff")
ta = read.table(fa, sep = "\t", as.is = T, header = F)

x = regexec("ID=([\\w\\.]+);Name=([\\w\\.]+);", ta$V9, perl = T)
begs = sapply(x, "[", 2)
lens = sapply(x, myfunc <- function(x) attr(x, "match.length")[2])
aids = substr(ta$V9, begs, begs + lens - 1)
begs = sapply(x, "[", 3)
lens = sapply(x, myfunc <- function(x) attr(x, "match.length")[3])
gids = substr(ta$V9, begs, begs + lens - 1)
ta2 = cbind.data.frame(ta[,c(1,4:7)], aid = aids, gid = gids, stringsAsFactors = F)
colnames(ta2)[1:5] = c("chr", "beg", "end", "score", "srd")

gb = group_by(ta2, aid)
ta3 = summarise(gb, chr = chr[1], beg = min(beg), end = max(end), score = mean(score), srd = srd[1], gid = gid[1])
colnames(ta3)[c(2:4,6)] = c("cid", "cbeg", "cend", "csrd")

tx = merge(ta3, ti2[,-2], by = 'cid', all.x = T)
fbegs = c(); fends = c(); fsrds = c()
for (i in 1:nrow(tx)) {
	if(is.na(tx$chr[i])) {
		fbegs = c(fbegs, NA); fends = c(fends, NA); fsrds = c(fsrds, NA)
		next
	}
	cbeg = tx$cbeg[i]; cend = tx$cend[i]; csrd = tx$csrd[i]
	beg = tx$beg[i]; end = tx$end[i]; srd = tx$srd[i]
	if(srd == "-") {
		fbeg = end - cend + 1
		fend = end - cbeg + 1
		fsrd = ifelse(csrd == "+", "-", "+")
	} else {
		fbeg = beg + cbeg - 1
		fend = beg + cend - 1
		fsrd = csrd
	}
	fbegs = c(fbegs, fbeg); fends = c(fends, fend); fsrds = c(fsrds, fsrd)
}

tx2 = cbind(tx[,c('gid','aid','score','cid','cbeg','cend','csrd','chr')], fbeg = fbegs, fend = fends, fsrd = fsrds, sid = tx$pid)
tx3 = tx2[order(tx2$gid),]
tx4 = tx3[tx3$sid %in% tp$V1,]

fo = file.path(dirw, "41.tsv")
write.table(tx4, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

