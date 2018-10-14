require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
source('Location.R')

dirw = file.path(Sys.getenv('misc2'), 'pbjoin')

##### prepare PBDT / PBBNDT join/break sites
alg = "PBDT"
alg = "PBBNDT"
fs = sprintf("%s/HM340.%s/raw.fix.fas.map", Sys.getenv("genome"), alg)
ts = read.table(fs, header = F, sep = "\t", as.is = T)
colnames(ts) = c("chr", "size", "chr2")
smap = ts$chr2; names(smap) = ts$chr

f01 = sprintf("%s/%s.txt", dirw, alg)
	t01 = read.table(f01, header = F, sep = "\t", as.is= T)
	colnames(t01) = c("nchr", "ochr", "obeg", "oend", "srd", "nbeg", "nend")

### 'joins'
x = table(t01$nchr)
tj1 = t01[t01$nchr %in% names(x)[x>1],]
tj2 = ddply(tj1, .(nchr), function(ds) {
	ds = ds[order(ds$nbeg),]
	ds = cbind(ds, type = "lr")
	ds$type = as.character(ds$type)
	ds$type[1] = 'r'
	ds$type[nrow(ds)] = 'l'
#	if(nrow(ds) > 2) {ds$type[2:(nrow(ds)-1)] = 'lr'}
	ds
})

chrs = c(); begs = c(); ends = c()
len = 10000
for (i in 1:nrow(tj2)) {
	if(tj2$type[i] %in% c('r', 'lr')) {
		chrs = c(chrs, tj2$nchr[i])
		begs = c(begs, max(tj2$nbeg[i]+1, tj2$nend[i]-len+1))
		ends = c(ends, tj2$nend[i])
	}
	if(tj2$type[i] %in% c('l', 'lr')) {
		chrs = c(chrs, tj2$nchr[i])
		begs = c(begs, tj2$nbeg[i]+1)
		ends = c(ends, min(tj2$nbeg[i]+len, tj2$nend[i]))
	}
}
tj = data.frame(chr = chrs, beg = begs, end = ends, stringsAsFactors = F)
grj = with(tj, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
grj = reduce(grj)
tj = data.frame(chr = as.character(seqnames(grj)), beg = start(grj), end = end(grj), stringsAsFactors = F)

### 'breaks'
x = table(t01$ochr)
tb1 = t01[t01$ochr %in% names(x)[x>1],]
tb1 = tb1[order(tb1$ochr, tb1$obeg),]
tb2 = ddply(tb1, .(ochr), function(ds) {
	ds = ds[order(ds$obeg),]
	ds = cbind(ds, type = "lr")
	ds$type = as.character(ds$type)
	ds$type[1] = 'r'
	ds$type[nrow(ds)] = 'l'
#	if(nrow(ds) > 2) {ds$type[2:(nrow(ds)-1)] = 'lr'}
	ds
})

chrs = c(); begs = c(); ends = c()
len = 10000
for (i in 1:nrow(tb2)) {
	if(tb2$type[i] %in% c('r', 'lr')) {
		chrs = c(chrs, tb2$nchr[i])
		begs = c(begs, max(tb2$nbeg[i]+1, tb2$nend[i]-len+1))
		ends = c(ends, tb2$nend[i])
	}
	if(tb2$type[i] %in% c('l', 'lr')) {
		chrs = c(chrs, tb2$nchr[i])
		begs = c(begs, tb2$nbeg[i]+1)
		ends = c(ends, min(tb2$nbeg[i]+len, tb2$nend[i]))
	}
}
tb = data.frame(chr = chrs, beg = begs, end = ends, stringsAsFactors = F)
grb = with(tb, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
grb = reduce(grb)
tb = data.frame(chr = as.character(seqnames(grb)), beg = start(grb), end = end(grb), stringsAsFactors = F)

tp = rbind(cbind(tj, type='join'), cbind(tb, type='break'))
tp = cbind(tp, len = tp$end-tp$beg+1)
tp$chr = smap[tp$chr]

fo = sprintf("%s/pbjoin/11.%s.tsv", Sys.getenv("misc2"), alg)
write.table(tp, fo, col.names = T, row.names = F, sep='\t', quote = F)

##### prepare PBBN / PBDTBN join sites
alg = "PBBN"
alg = "PBDTBN"

fs = sprintf("%s/HM340.%s/raw.fix.fas.map", Sys.getenv("genome"), alg)
ts = read.table(fs, header = F, sep = "\t", as.is = T)
colnames(ts) = c("chr", "size", "chr2")
smap = ts$chr2; names(smap) = ts$chr

f01 = sprintf("%s/%s_join.txt", dirw, alg)
	t01 = read.table(f01, header = T, sep = "\t", as.is= T)
	t01 = t01[t01$Compnt_Type=='W', c(1,6:9,2,3)]
	colnames(t01) = c("nchr", "ochr", "obeg", "oend", "srd", "nbeg", "nend")
	t01$nchr = smap[t01$nchr]

### 'joins'
x = table(t01$nchr)
tj1 = t01[t01$nchr %in% names(x)[x>1],]
tj2 = ddply(tj1, .(nchr), function(ds) {
	ds = ds[order(ds$nbeg),]
	ds = cbind(ds, type = "lr")
	ds$type = as.character(ds$type)
	ds$type[1] = 'r'
	ds$type[nrow(ds)] = 'l'
#	if(nrow(ds) > 2) {ds$type[2:(nrow(ds)-1)] = 'lr'}
	ds
})

chrs = c(); begs = c(); ends = c()
len = 10000
for (i in 1:nrow(tj2)) {
	if(tj2$type[i] %in% c('r', 'lr')) {
		chrs = c(chrs, tj2$nchr[i])
		begs = c(begs, max(tj2$nbeg[i]+1, tj2$nend[i]-len+1))
		ends = c(ends, tj2$nend[i])
	}
	if(tj2$type[i] %in% c('l', 'lr')) {
		chrs = c(chrs, tj2$nchr[i])
		begs = c(begs, tj2$nbeg[i]+1)
		ends = c(ends, min(tj2$nbeg[i]+len, tj2$nend[i]))
	}
}
tj = data.frame(chr = chrs, beg = begs, end = ends, stringsAsFactors = F)
grj = with(tj, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
grj = reduce(grj)
tj = data.frame(chr = as.character(seqnames(grj)), beg = start(grj), end = end(grj), stringsAsFactors = F)

tp = cbind(tj, type='join')
tp = cbind(tp, len = tp$end-tp$beg+1)
fo = sprintf("%s/pbjoin/11.%s.tsv", Sys.getenv("misc2"), alg)
write.table(tp, fo, col.names = T, row.names = F, sep='\t', quote = F)

#### generate background comparison by sampling genome windows
algs = c("PBBN", "PBDT", "PBBNDT", "PBDTBN")
for (alg in algs) {
fs = sprintf("%s/HM340.%s/raw.fix.fas.map", Sys.getenv("genome"), alg)
ts = read.table(fs, header = F, sep = "\t", as.is = T)
colnames(ts) = c("chr", "size", "chr2")

seqlens = ts$size
names(seqlens) = ts$chr2
grw = unlist(tileGenome(seqlens, tilewidth = 10000))
grw = grw[width(grw) > 5000]
tw = data.frame(chr = as.character(seqnames(grw)), beg = start(grw), end = end(grw), stringsAsFactors = F)
tw = cbind(tw, len = tw$end - tw$beg + 1)
fo = sprintf("%s/pbjoin/15.bg.%s.tsv", Sys.getenv("misc2"), alg)
write.table(tw, fo, col.names = T, row.names = F, sep='\t', quote = F)
}

#### look for BN bias (PBBNDT)
alg = 'PBBNDT'

fp = sprintf("%s/pbjoin/11.%s.tsv", Sys.getenv("misc2"), alg)
tp = read.table(fp, header = T, sep = "\t", as.is= T)
grp = with(tp, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

fw = sprintf("%s/pbjoin/15.bg.%s.tsv", Sys.getenv("misc2"), alg)
tw = read.table(fw, header = T, sep = "\t", as.is= T)
grw = with(tw, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

res = "BbvCI"
f02 = sprintf("%s/02.restriction/HM340.%s.%s.txt", dirw, alg, res)
t02 = read.table(f02, header = T, sep = "\t", as.is= T)
colnames(t02) = c("chr", "motif", "beg", "end", "mm", "srd", "str")
grr1 = with(t02, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

res = "BspQI"
f02 = sprintf("%s/02.restriction/HM340.%s.%s.txt", dirw, alg, res)
t02 = read.table(f02, header = T, sep = "\t", as.is= T)
colnames(t02) = c("chr", "motif", "beg", "end", "mm", "srd", "str")
grr2 = with(t02, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

cnt1 = intersect_count(grp, grr1)
cnt2 = intersect_count(grp, grr2)
tp = cbind(tp, den_bbv = cnt1/tp$len*1000, den_bsp = cnt2/tp$len*1000)

cnt1 = intersect_count(grw, grr1)
cnt2 = intersect_count(grw, grr2)
tw = cbind(tw, den_bbv = cnt1/tw$len*1000, den_bsp = cnt2/tw$len*1000)

to1 = rbind(tp[,c('type','den_bbv')], data.frame(type='random', den_bbv=tw$den_bbv, stringsAsFactors = F))
to1 = cbind(to1, enz = 'BbvCI')
colnames(to1)[2] = 'den'
to2 = rbind(tp[,c('type','den_bsp')], data.frame(type='random', den_bsp=tw$den_bsp, stringsAsFactors = F))
to2 = cbind(to2, enz = 'BspQI')
colnames(to2)[2] = 'den'
to = rbind(to1, to2)
p1 = ggplot(to) +
  geom_density(aes(den, fill = type), alpha = 0.5, adjust = 4) + 
  scale_x_continuous(name = '# restriction enzyme cut sites (per 1kb)', limits=c(0,1.5)) +
  scale_y_continuous(name = 'Density') +
  scale_fill_brewer(palette = 'Accent') +
  facet_grid(. ~ enz) +
  theme_bw() +
  ggtitle(alg) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(legend.position = 'top', legend.background = element_blank(), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(1, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_blank(), legend.text = element_text(size = 9, angle = 0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = sprintf("%s/31.%s.pdf", dirw, alg)
ggsave(p1, filename = fp, width = 7, height = 5)

y1 = to$den[to$type=='join' & to$enz=='BbvCI']
y2 = to$den[to$type=='break' & to$enz=='BbvCI']
y = to$den[to$type=='random' & to$enz=='BbvCI']
t.test(y1, y, alternative = 'greater')
t.test(y2, y, alternative = 'greater')
t.test(c(y1,y2), y, alternative = 'greater')

y1 = to$den[to$type=='join' & to$enz=='BspQI']
y2 = to$den[to$type=='break' & to$enz=='BspQI']
y = to$den[to$type=='random' & to$enz=='BspQI']
t.test(y1, y, alternative = 'greater')
t.test(y2, y, alternative = 'greater')
t.test(c(y1,y2), y, alternative = 'greater')

#### look for DT bias (PBDTBN)
alg = "PBDTBN"

fp = sprintf("%s/pbjoin/11.%s.tsv", Sys.getenv("misc2"), alg)
tp = read.table(fp, header = T, sep = "\t", as.is= T)
grp = with(tp, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

fw = sprintf("%s/pbjoin/15.bg.%s.tsv", Sys.getenv("misc2"), alg)
tw = read.table(fw, header = T, sep = "\t", as.is= T)
grw = with(tw, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

res = "DpnII"
f03 = sprintf("%s/02.restriction/HM340.%s.%s.txt", dirw, alg, res)
t03 = read.table(f03, header = T, sep = "\t", as.is= T)
colnames(t03) = c("chr", "motif", "beg", "end", "mm", "srd", "str")
grr3 = with(t03, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

cnt = intersect_count(grp, grr3)
tp = cbind(tp, den_dpn = cnt/tp$len*1000)

cnt = intersect_count(grw, grr3)
tw = cbind(tw, den_dpn = cnt/tw$len*1000)

to = rbind(tp[,c('type','den_dpn')], data.frame(type='random', den_dpn=tw$den_dpn, stringsAsFactors = F))
to = cbind(to, enz = 'DpnII')
colnames(to)[2] = 'den'
p1 = ggplot(to) +
  geom_density(aes(den, fill = type), alpha = 0.5, adjust = 2) + 
  scale_x_continuous(name = '# restriction enzyme cut sites (per 1kb)') +
  scale_y_continuous(name = 'Density') +
  scale_fill_brewer(palette = 'Accent') +
  facet_grid(. ~ enz) +
  theme_bw() +
  ggtitle(alg) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(legend.position = 'top', legend.background = element_blank(), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(1, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_blank(), legend.text = element_text(size = 9, angle = 0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = sprintf("%s/31.%s.pdf", dirw, alg)
ggsave(p1, filename = fp, width = 5, height = 5)

y1 = to$den[to$type=='join' & to$enz=='DpnII']
y = to$den[to$type=='random' & to$enz=='DpnII']
t.test(y1, y, alternative = 'greater')



#### look for repeat enrichment
#alg = "PBBN"
#alg = "PBDT"
alg = "PBBNDT"
#alg = "PBDTBN"

fp = sprintf("%s/pbjoin/11.%s.tsv", Sys.getenv("misc2"), alg)
tp = read.table(fp, header = T, sep = "\t", as.is= T)
grp = with(tp, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

fw = sprintf("%s/pbjoin/15.bg.%s.tsv", Sys.getenv("misc2"), alg)
tw = read.table(fw, header = T, sep = "\t", as.is= T)
grw = with(tw, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

fm = sprintf("%s/HM340.%s/12.rm.bed", Sys.getenv('genome'), alg)
tm = read.table(fm, header = F, as.is = T, sep = "\t")
grm = with(tm, GRanges(seqnames = V1, ranges = IRanges(V2+1, end = V3)))

cnt = intersect_basepair(grp, grm)
tp = cbind(tp, prop = cnt/tp$len)

cnt = intersect_basepair(grw, grm)
tw = cbind(tw, prop = cnt/tw$len)

to = rbind(tp[,c('type','prop')], data.frame(type='random', prop=tw$prop, stringsAsFactors = F))
p1 = ggplot(to) +
  geom_density(aes(prop, fill = type), alpha = 0.5) + 
  scale_x_continuous(name = 'Repeat Content', limits = c(0,1)) +
  scale_y_continuous(name = 'Density') +
  scale_fill_brewer(palette = 'Accent') +
  theme_bw() +
  ggtitle(alg) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(legend.position = 'top', legend.background = element_blank(), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(1, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_blank(), legend.text = element_text(size = 9, angle = 0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = sprintf("%s/51.repeat.%s.pdf", dirw, alg)
ggsave(p1, filename = fp, width = 5, height = 5)

y1 = to$prop[to$type=='join']
y2 = to$prop[to$type=='break']
y = to$prop[to$type=='random']
t.test(y1, y, alternative = 'greater')
t.test(y2, y, alternative = 'greater')
t.test(c(y1,y2), y, alternative = 'greater')

### prepare PBBN/PBDT join/break coordinates
dirw = file.path(Sys.getenv('misc1'), 'hm340.ms')
fh = file.path(dirw, "00.shared.tsv")
th = read.table(fh, header = F, sep = "\t", as.is = T)

alg = "PBBN"
#alg = "PBDT"
fs = sprintf("%s/HM340.%s/raw.fix.fas.map", Sys.getenv("genome"), alg)
ts = read.table(fs, header = F, sep = "\t", as.is = T)
colnames(ts) = c("chr", "size", "chr2")
smap = ts$chr2; names(smap) = ts$chr

f01 = sprintf("%s/01.%s.tsv", dirw, tolower(alg))
t01 = read.table(f01, header = F, sep = "\t", as.is= T)
colnames(t01) = c("chr", "beg", "end")
t01 = cbind.data.frame(t01, shared = '', stringsAsFactors = F)
if(alg == "PBBN") {
	t01$shared[t01$chr %in% th$V1] = 1
} else {
	t01$shared[t01$chr %in% th$V2] = 1
}
stopifnot(t01$chr %in% ts$chr)
t01$chr = smap[t01$chr]

if(alg == "PBDT") t01$beg = t01$beg + 1
fo = sprintf("%s/03.%s.tsv", dirw, tolower(alg))
write.table(t01, fo, col.names = F, row.names = F, sep='\t', quote = F)
