require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)

org = 'W22'
org = 'PH207'
dirw = sprintf("/home/springer/zhoux379/data/genome/%s/synteny", org)

# reference gene ordering
fg = '/home/springer/zhoux379/data/genome/B73/52.rep.gtb'
tg = read.table(fg, header = T, sep = "\t", as.is = T)[,c(2:6)]
tg = tg[order(tg$chr, (tg$beg + tg$end)/2),]
tg = cbind(tg, num = NA)
for (chr in unique(tg$chr)) {
	tg$num[tg$chr == chr] = 1:sum(tg$chr == chr)
}
tg = cbind(tg, pos = as.integer((tg$beg+tg$end)/2))
tgt = tg[,c(1,2,7,5,6)]
colnames(tgt) = c("tid", "tchr", "tpos", "tsrd", "tnum")

# query gene ordering
fg = sprintf("%s/../52.rep.gtb", dirw)
tg = read.table(fg, header = T, sep = "\t", as.is = T)[,c(2:6)]
tg = tg[order(tg$chr, (tg$beg + tg$end)/2),]
tg = cbind(tg, num = NA)
for (chr in unique(tg$chr)) {
	tg$num[tg$chr == chr] = 1:sum(tg$chr == chr)
}
tg = cbind(tg, pos = as.integer((tg$beg+tg$end)/2))
tgq = tg[,c(1,2,7,5,6)]
colnames(tgq) = c("qid", "qchr", "qpos", "qsrd", "qnum")

# read blastp scores
fi = sprintf("%s/%s_B73.tsv", dirw, org)
ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1,3,7)]
colnames(ti) = c("qid", "tid", "evalue")
ti$evalue[ti$evalue == 0] = min(ti$evalue[ti$evalue != 0])
min(ti$evalue)
ti = cbind(ti, logE = -log(ti$evalue))

# read gene mapping
fm = sprintf("%s/mapping.tsv", dirw)
tm = read.table(fm, header = F, sep = "\t", as.is = T)[,c(4,8)]
colnames(tm) = c("tid", "qid")

tm1 = merge(tm, tgq[,c(1,2,5)], by = c("qid"))
nrow(tm1) - nrow(tm)

tm2 = merge(tm1, tgt[,c(1,2,5)], by = c("tid"))
nrow(tm2) - nrow(tm1)

tm3 = merge(tm2, ti[,-3], by = c("qid", "tid"))
grp = group_by(tm3, qid, tid)
tm4 = summarise(grp, qchr = qchr[1], tchr = tchr[1], qnum = qnum[1], tnum = tnum[1], logE = max(logE))
nrow(tm4) - nrow(tm2)

# output for dagchainer

# output for quota-align
to = tm4[,c('qchr','qnum','tchr','qnum','logE')]
fo = file.path(dirw, "10.qa")
write.table(to, fo, sep = "\t", row.names = F, col.names = F, quote = F)
# python2 $src/git/quota-alignment/quota_align.py --format=raw --merge --Dm=20 --min_size=5 --quota=1:1 10.qa 



## plot
chrmap = c()
for (i in 1:10) {
    chrmap[sprintf("%s", i)] = i
    chrmap[sprintf("chr%d", i)] = i
    chrmap[sprintf("chr%02d", i)] = i
}

# read in ref/query chromosome sizes
fz = sprintf('/home/springer/zhoux379/data/genome/%s/15.sizes', "B73")
tz = read.table(fz, sep="\t", header=F, as.is=T)
colnames(tz)=c("chr","size")
tz = tz[tz$chr %in% names(chrmap),]
tz = cbind(tz, chrnum = chrmap[tz$chr])
tz = tz[order(tz$chrnum),]
offsets = c(0, cumsum(tz$size + 20000000)[-nrow(tz)])
tz = cbind(tz, offset = offsets)
tz = within(tz, {beg = 1+offset; end = size+offset; pos =(beg+end)/2})
tzt = tz
colnames(tzt)[c(1,4)] = c("tchr", 'toffset')

fz = sprintf('/home/springer/zhoux379/data/genome/%s/15.sizes', org)
tz = read.table(fz, sep="\t", header=F, as.is=T)
colnames(tz)=c("chr","size")
tz = tz[tz$chr %in% names(chrmap),]
tz = cbind(tz, chrnum = chrmap[tz$chr])
tz = tz[order(tz$chrnum),]
offsets = c(0, cumsum(tz$size + 20000000)[-nrow(tz)])
tz = cbind(tz, offset = offsets)
tz = within(tz, {beg = 1+offset; end = size+offset; pos =(beg+end)/2})
tzq = tz
colnames(tzq)[c(1,4)] = c("qchr", 'qoffset')

# dotplot
fq = sprintf("%s/10.qa.filtered", dirw)
tq = read.table(fq, header = F, sep = "\t", as.is = T)
colnames(tq) = c("qchr", "qnum", "tchr", "tnum", "logE")

tq2 = merge(tq, tgt, by = c('tchr','tnum'))
nrow(tq2) - nrow(tq)
tq3 = merge(tq2, tgq, by = c('qchr','qnum'))
nrow(tq3) - nrow(tq2)

tp = tq3

tp = tp[tp$tchr %in% names(chrmap) & tp$qchr %in% names(chrmap),]
nrow(tp) - nrow(tq3)

tp1 = merge(tp, tzt[,c(1,4)], by = 'tchr')
tp2 = merge(tp1, tzq[,c(1,4)], by = 'qchr')
tp3 = within(tp2, {tpos = tpos + toffset; qpos = qpos + qoffset; rm(qoffset, toffset)})

p1 = ggplot(tp3) +
  geom_point(aes(x = tpos, y = qpos), size=0.03) +
  geom_vline(xintercept = tzt$beg, alpha=0.3) +
  geom_vline(xintercept = tzt$end, alpha=0.3) +
  geom_hline(yintercept = tzq$beg, alpha=0.3) +
  geom_hline(yintercept = tzq$end, alpha=0.3) +
  scale_x_continuous(name=sprintf("%s","B73"), breaks=tzt$pos, labels=tzt$tchr, expand=c(0,0)) +
  scale_y_continuous(name=sprintf("%s",org), breaks=tzq$pos, labels=tzq$qchr, expand=c(0,0)) +
  #scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(panel.grid = element_blank()) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 9), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, color = "black", angle = 0, hjust = 1))
fp = sprintf("%s/15.filtered.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)

