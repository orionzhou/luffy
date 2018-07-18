require(Hmisc)
require(tidyr)
require(RColorBrewer)

dirw = "/home/springer/zhoux379/data/misc1/li2013"
chrs = sprintf("%d", 1:10)


fi = file.path(dirw, "eQTL_Li2013.tsv")
ti = read.table(fi, sep="\t", header=T, as.is=T)

## check how cis/trans is defined
ti2 = cbind(ti, samchr = ti$gchr.p == ti$qchr.p)
ti3 = ti2[ti2$type == 'cis',]
describe(abs(ti3$gpos.p - ti3$qpos.p))
ti4 = ti3[abs(ti3$gpos.p - ti3$qpos.p) > 10000000,]
describe(abs(ti4$gpos.g - ti4$qpos.g))

tig = within(ti, {
	type2 = ifelse(gchr.p == qchr.p,
	ifelse(abs(gpos.p-qpos.p) <= 10000000, 'cis', 
		ifelse(abs(gpos.g-qpos.g) <= 10, 'cis', 'trans')),
	'trans')
})
identical(tig$type2, ti$type)

x = tig[tig$type != tig$type2,]
x

## convert eQTL coordinates from v2 to v4
tb = data.frame(chr = ti$qchr.p, beg = ti$qpos.p-1, end = ti$qpos.p, stringsAsFactors=F)
tb = unique(tb[order(tb$chr, tb$beg),])
fo = file.path(dirw, '03.v2.bed')
write.table(tb, fo, sep = "\t", row.names = F, col.names = F, quote = F)

#CrossMap.py bed $genome/Zmays_v4/chain/AGPv2_to_AGPv4.chain.gz 03.v2.bed 04.v4.bed

fm = file.path(dirw, "04.v4.bed")
tm = read.table(fm, sep="\t", header=F, as.is=T)
colnames(tm) = c("nchr", "nbeg", "nend")
fu = file.path(dirw, "04.v4.bed.unmap")
tu = read.table(fu, sep="\t", header=F, as.is=T)
colnames(tu) = c("chr", "beg", "end")
stopifnot(nrow(tb) == nrow(tm) + nrow(tu))

tu = cbind(tu, unmap = 1)
tb2 = merge(tb, tu, by.x = c('chr','beg','end'), all.x = T)
tb3 = tb2[is.na(tb2$unmap),]
stopifnot(nrow(tb3) == nrow(tm))
tb4 = cbind(tb3[,-4], tm)
tb5 = tb4[,c('chr','end','nchr','nend')]
colnames(tb5) = c("chr","pos","nchr",'npos')

## read gene mapping
fm = '/home/springer/zhoux379/data/genome/Zmays_v4/gene_mapping/maize.v3TOv4.geneIDhistory.txt'
tm = read.table(fm, sep="\t", header=F, as.is=T)
colnames(tm) = c("gid", "ngid", "note", "method", "type")
gids = unique(ti$gid)
length(gids)
sum(gids %in% tm$gid)
tm2 = tm[tm$gid %in% gids,]
table(tm2$type)
tm3 = tm2[tm2$type == '1-to-1',]
stopifnot(nrow(tm3)==length(unique(tm3$gid)))

## conversion
ti2 = merge(ti, tm3[,1:2], by = 'gid')
ti3 = merge(ti2, tb5, by.x = c("qchr.p","qpos.p"), by.y = c("chr","pos"))
ti4 = ti3[order(ti3$ngid), c(4,18,9:13,19:20,14:17)]

fg = '/home/springer/zhoux379/data/genome/Zmays_v4/51.gtb'
tg = read.table(fg, sep="\t", header=T, as.is=T)
grp = dplyr::group_by(tg, par)
tg2 = dplyr::summarise(grp, gchr = chr[1], gbeg = min(beg), gend = max(end))

ti5 = merge(ti4, tg2, by.x = 'ngid', by.y = 'par')
ti5 = ti5[,c(2,1,14:16,3:13)]
colnames(ti5)[c(2,11,12)] = c("gid", "qchr", "qpos")

fo = file.path(dirw, '10.eQTL.v4.tsv')
write.table(ti5, fo, sep = "\t", row.names = F, col.names = T, quote = F)


### visualize
chrs = sprintf("%d", 1:10)

fi = file.path(dirw, "10.eQTL.v4.tsv")
ti = read.table(fi, sep="\t", header=T, as.is=T)
ti = ti[ti$gchr %in% chrs & ti$qchr %in% chrs,]
ti$gchr = as.integer(ti$gchr)
ti$qchr = as.integer(ti$qchr)


fz = '/home/springer/zhoux379/data/genome/Zmays_v4/15.sizes'
tz = read.table(fz, sep="\t", header=F, as.is=T)
colnames(tz)=c("chr","size")
tz = tz[tz$chr %in% chrs,]
tz$chr = as.integer(tz$chr)
offsets = c(0, cumsum(tz$size + 20000000)[-nrow(tz)])
tz = cbind(tz, offset = offsets)
tz = within(tz, {beg = 1+offset; end = size+offset; pos =(beg+end)/2})

ti1 = merge(ti, tz[,c(1,3)], by.x = 'gchr', by.y = 'chr')
ti1 = within(ti1, {gpos = gbeg + offset; rm(offset)})
ti2 = merge(ti1, tz[,c(1,3)], by.x = 'qchr', by.y = 'chr')
ti2 = within(ti2, {qpos = qpos + offset; rm(offset)})

tp = ti2

p1 = ggplot(tp) +
  geom_point(aes(x = qpos, y = gpos, color =type), size=0.2) +
  geom_vline(xintercept = tz$beg, alpha=0.3) +
  geom_vline(xintercept = tz$end, alpha=0.3) +
  geom_hline(yintercept = tz$beg, alpha=0.3) +
  geom_hline(yintercept = tz$end, alpha=0.3) +
  scale_x_continuous(name='eQTL position', breaks=tz$pos, labels=tz$chr, expand=c(0,0)) +
  scale_y_continuous(name='eGene position', breaks=tz$pos, labels=tz$chr, expand=c(0,0)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(panel.grid = element_blank()) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 9), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, color = "black", angle = 0, hjust = 1))
fp = sprintf("%s/11.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)


### identify hotspots - sliding window approach
fi = file.path(dirw, "10.eQTL.v4.tsv")
ti = read.table(fi, sep="\t", header=T, as.is=T)
ti = within(ti, {type2 = ifelse(gchr==qchr & abs((gbeg+gend)/2-qpos) < 10000000, 'cis', 'trans')})

tq = data.frame(chr=ti$qchr, beg=ti$qpos-1, end=ti$qpos, id=ti$id, gid=ti$gid, type=ti$type2)
tq = tq[order(tq$chr,tq$beg),]

fo = file.path(dirw, '22.eqtl.bed')
write.table(tq, fo, sep = "\t", row.names = F, col.names = F, quote = F)

# cd $misc1/li2013
# winslide.pl -i $genome/Zmays_v4/15.sizes -o 21.win.bed -step 1000000 -size 1000000
# intersectBed -wao -a 21.win.bed -b 22.eqtl.bed > 23.ovlp.bed
# intersectBed -wao -a 21.win.bed -b $genome/Zmays_v4/v37/gene.bed > 21.win.gene.bed

fv = file.path(dirw, "23.ovlp.bed")
tv = read.table(fv, sep="\t", header=F, as.is=T)[,1:9]
colnames(tv) = c("chr","beg","end",'qchr','qbeg','qend','id','gid','type')

tv = tv[tv$chr %in% chrs,]
tv$chr = as.integer(tv$chr)

grp = dplyr::group_by(tv, chr, beg, end)
tv2 = dplyr::summarise(grp, cis = sum(type=='cis'), trans=sum(type=='trans'), nsnp = length(unique(qend[type=='trans'])))
tv3 = merge(tv2, tz[,c('chr','offset')], by = 'chr')
tv3 = within(tv3, {pos = (beg+end)/2+offset+1})
tv4 = gather(tv3[,c("chr","pos",'cis','trans')], type, ngene, -chr, -pos)
tv4$type = factor(tv4$type, levels = c("trans", "cis"))

p1 = ggplot(tv4) +
  geom_bar(aes(x = pos, y = ngene, fill = type), stat='identity') +
  geom_segment(data=tz, aes(x=beg,xend=end,y=0,yend=0), size=1.5) +
  #geom_vline(xintercept = tz$beg, alpha=0.3) +
  #geom_vline(xintercept = tz$end, alpha=0.3) +
  scale_x_continuous(name='eQTL position', breaks=tz$pos, labels=tz$chr, expand=c(0,0)) +
  scale_y_continuous(name='# of eGenes', expand=c(0,0)) +
  scale_fill_manual(values = brewer.pal(3, "Set1")[2:1]) +
  facet_wrap(~type, nrow=2) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(panel.grid = element_blank()) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'none', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, color = "black", angle = 0, hjust = 1))
fp = sprintf("%s/23.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 6)



fg = file.path(dirw, '21.win.gene.bed')
tg = read.table(fg, sep="\t", header=F, as.is=T)
colnames(tg) = c("chr","beg","end","gchr","gbeg","gend","gid",'bp')
tg = tg[tg$chr %in% chrs,]
tg$chr = as.integer(tg$chr)

grp = dplyr::group_by(tg, chr, beg, end)
tg2 = dplyr::summarise(grp, ngene = sum(gid != '.'))
stopifnot(nrow(tg2)==nrow(tv2))

tv5 = cbind.data.frame(tv2, ngene = tg2$ngene)
tv5 = within(tv5, {den = trans/(1*ngene)})
tv5$den[is.na(tv5$den)] = 0

tv6 = tv5[tv5$den >= 1.25 & tv5$trans >= 10,]

fo = file.path(dirw, "25.hs.raw.bed")
write.table(tv6[,1:3], fo, sep = "\t", row.names = F, col.names = F, quote = F)

# mergeBed -i 25.hs.raw.bed > 26.hs.merged.bed
# intersectBed -wao -a 26.hs.merged.bed -b 25.hs.raw.bed > 27.hs.bed

fm = file.path(dirw, '27.hs.bed')
tm = read.table(fm, sep="\t", header=F, as.is=T)[,1:6]
colnames(tm) = c("mchr","mbeg","mend","chr","beg","end")
tms = unique(tm[,1:3])
tms = tms[order(tms$mchr, tms$mbeg),]
tms = cbind(tms, mid = sprintf("hs%02d", 1:nrow(tms)))
tm = merge(tms, tm, by=c("mchr","mbeg","mend"))


tm2 = merge(tm, tv, by=c("chr","beg","end"))[,-c(1:3)]
colnames(tm2)[1:3] = c("chr","beg","end")
tm2 = tm2[order(tm2$chr,tm2$beg),]
fo = file.path(dirw, "29.hs.tsv")
write.table(tm2, fo, sep = "\t", row.names = F, col.names = T, quote = F)

tm3 = tm2[tm2$type=='trans',]
tm4 = data.frame(mid = sprintf("trans-eQTL.%s", tm3$mid), gid = tm3$gid, funcat='', opt='trans-eQTL', stringsAsFactors = F)

fo = file.path(dirw, "29.hs.module.tsv")
write.table(tm4, fo, sep = "\t", row.names = F, col.names = T, na = '', quote = F)

# generate stats
fh = file.path(dirw, "29.hs.tsv")
th = read.table(fh, sep="\t", header=T, as.is=T)
ths = ddply(th, .(mid), summarise, ntrans = sum(type=='trans'))



### identify trans-eQTL hotspots - single SNP approach
fi = file.path(dirw, "10.eQTL.v4.tsv")
ti = read.table(fi, sep="\t", header=T, as.is=T)
ti = within(ti, {type2 = ifelse(gchr==qchr & abs((gbeg+gend)/2-qpos) < 10000000, 'cis', 'trans')})

tq = unique(data.frame(chr=ti$qchr, beg=ti$qpos-1, end=ti$qpos))
tq = tq[order(tq$chr,tq$beg),]

fo = file.path(dirw, '41.eqtl.bed')
write.table(tq, fo, sep = "\t", row.names = F, col.names = F, quote = F)

# bcftools view -H -R 41.eqtl.bed $misc2/mo17vnt/53.vnt.final/02.vcf.gz > 42.vnt.vcf

ti2 = ti[ti$type2 == 'trans',]
ti2 = cbind(ti, direction = sign(ti$ADD))
grp = dplyr::group_by(ti2, qchr, qpos, direction)
ti3 = dplyr::summarise(grp, ntgt = n())
sum(ti3$ntgt >= 20)
ti4 = ti3[ti3$ntgt >= 20,]

ti4 = ti4[ti4$qchr %in% chrs,]
ti4$qchr = as.integer(ti4$qchr)
ti4 = ti4[order(ti4$qchr, ti4$qpos),]
ti4a = ti4[ti4$direction == -1,]
ti4a = cbind(mid = sprintf("trans-eQTL-B.hs%02d", 1:nrow(ti4a)), ti4a)
ti4b = ti4[ti4$direction == 1,]
ti4b = cbind(mid = sprintf("trans-eQTL-M.hs%02d", 1:nrow(ti4b)), ti4b)
ti4 = rbind(ti4a, ti4b)


ti5 = merge(ti2, ti4[,-5], by=c("qchr","qpos",'direction'))
ti5$qchr = as.integer(ti5$qchr)
ti5 = ti5[order(ti5$mid, ti5$gid),]

fo = file.path(dirw, "49.hs.tsv")
write.table(ti5, fo, sep = "\t", row.names = F, col.names = T, na = '', quote = F)

tm = data.frame(mid = ti5$mid, gid = ti5$gid, funcat='', opt='trans-eQTL', stringsAsFactors = F)
fo = file.path(dirw, "49.hs.module.tsv")
write.table(tm, fo, sep = "\t", row.names = F, col.names = T, na = '', quote = F)

# generate stats
fh = file.path(dirw, "49.hs.module.tsv")
th = read.table(fh, sep="\t", header=T, as.is=T)
ths = ddply(th, .(mid), summarise, size=length(mid))
describe(ths$size)
