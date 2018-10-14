require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(dplyr)
source("comp.fun.R")
source("Location.R")

dirw = file.path(Sys.getenv("misc3"), "comp.vnt")
diro = file.path(Sys.getenv("misc3"), "comp.genefam")

tg = read.table(tcfg$gene, sep = "\t", header = F, as.is = T)
colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
tg = tg[tg$type == 'cds',]
grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

##### calculate ThetaPi for each gene
fa = file.path(dirw, "81.cvg.tbl")
ta = read.table(fa, header = F, sep = "\t", as.is = T)
gra = with(ta, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

fv = file.path(dirw, '25.stat.tbl')
tv = read.table(fv, header = T, sep = "\t", as.is = T)
gsnp = with(tv, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))

effmap = c(
  'downstream_gene_variant' = 'intergenic',
  'upstream_gene_variant' = 'intergenic',
  'intergenic_region' = 'intergenic',
  'missense_variant' = 'replacement',
  'synonymous_variant' = 'synonymous',
  'stop_retained_variant' = 'synonymous',
  'initiator_codon_variant' = 'synonymous',
  'intron_variant' = 'intron',
  'splice_region_variant' = 'intron',
  '5_prime_UTR_variant' = 'utr',
  '3_prime_UTR_variant' = 'utr',
  '5_prime_UTR_premature_start_codon_gain_variant' = 'utr',
  'stop_gained' = 'large_effect',
  'splice_acceptor_variant' = 'large_effect',
  'splice_donor_variant' = 'large_effect',
  'stop_lost' = 'large_effect',
  'start_lost' = 'large_effect'
)
ff = file.path(dirw, "23.snp.anno.tbl")
tf = read.table(ff, header = F, sep = "\t", as.is = T)
stopifnot(nrow(tv) == nrow(tf))
stopifnot(sum(!tf$V3 %in% names(effmap)) == 0)
tv = cbind(tv, eff = as.character(effmap[tf$V3]))


nds = intersect_score(grg, gsnp)
bps = intersect_basepair(grg, gra)

tg2 = cbind(tg, nd = nds, bp = bps)
gb = group_by(tg2, id)
dg = dplyr::summarise(gb, fam = fam[1], len = sum(end-beg+1), lenc = sum(bp), nd = sum(nd))

dg = dg[dg$lenc / dg$len >= 0.8,]
genes_covered = dg$id
to = cbind(dg[,c('id','fam')], pi = dg$nd / dg$lenc)

fo = file.path(diro, "03.pi.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

tp = ddply(to, .(fam), summarise, cnt = length(fam), q5=quantile(pi,0.05), q25=quantile(pi,0.25), q50=quantile(pi,0.5), q75=quantile(pi,0.75), q95=quantile(pi,0.95))
tp = tp[order(tp$q50, decreasing=T),]
fams = tp$fam
tp$fam = factor(tp$fam, levels = fams)

p1 = ggplot(tp) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', position = 'dodge', width = 0.7) +
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01), labels = sprintf("%s | %5d", tp$fam, tp$cnt)) +
  scale_y_continuous(name = 'ThetaPi', expand = c(0, 0), limits = c(0, 0.02)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0,'lines')) +
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = file.path(diro, "03.pi.pdf")
ggsave(p1, filename = fp, width = 5, height = 8)

### characterize large-effect SNPs by gene fam
lefmap = c(
  'stop_gained' = 'premature stop',
  'splice_acceptor_variant' = 'splice site',
  'splice_donor_variant' = 'splice site',
  'stop_lost' = 'stop codon lost',
  'start_lost' = 'start codon lost'
)
idxs = tf$V3 %in% names(lefmap)
tvl = cbind(tv[idxs, c('chr','pos')], eff = as.character(lefmap[tf$V3[idxs]]))

dz = data.frame()
for (lef in lefmap) {
  tvls = tvl[tvl$eff == lef,]
  gvl = with(tvls, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))
  cnts = intersect_count(grg, gvl)
  dz = rbind(dz, cbind(tg, eff = lef, cnt = cnts))
}
gb = group_by(dz, id, eff)
dz2 = dplyr::summarize(gb, fam = fam[1], cnt = sum(cnt))
to = dz2[dz2$id %in% genes_covered,]


fo = file.path(diro, "07.largeeff.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

tp = ddply(to, .(fam, eff), summarise, prop = sum(cnt > 0) / length(cnt))
tp2 = ddply(tp, .(fam), summarise, prop = sum(prop))
fams = tp2$fam[order(tp2$prop, decreasing=T)]
tp$fam = factor(tp$fam, levels = fams)

p2 = ggplot(tp) +
  geom_bar(aes(x = fam, y = prop, fill = eff),
    stat = 'identity', position = 'stack', width = 0.7) + 
  coord_flip() +
  scale_x_discrete(name = '', breaks = fams, labels = fams, expand = c(0.01, 0.01)) +
  scale_y_continuous(name = 'Proportion w. large-effect changes', expand = c(0, 0), limits = c(0, 1)) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  theme(legend.position = c(0.7, 0.8), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.7, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_blank(), legend.text = element_text(size = 8, angle = 0)) +
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = file.path(diro, "07.largeeff.pdf")
ggsave(p2, filename = fp, width = 5, height = 8)

### Gene-Family Theta-Pi using reference mapping-based approach
dirw = file.path(Sys.getenv("misc3"), "hapmap/12_ncgr")
diro = file.path(Sys.getenv("misc3"), "comp.stat")

fa = file.path(dirw, "81.cvg.tbl")
ta = read.table(fa, header = F, sep = "\t", as.is = T)
gra = with(ta, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

fv = file.path(dirw, '45.stat.tbl')
tv = read.table(fv, header = T, sep = "\t", as.is = T)
gsnp = with(tv, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))

nds = intersect_score(grg, gsnp)
bps = intersect_basepair(grg, gra)

tg2 = cbind(tg, nd = nds, bp = bps)
gb = group_by(tg2, id)
dg = dplyr::summarise(gb, fam = fam[1], len = sum(end-beg+1), lenc = sum(bp), nd = sum(nd))

dg = dg[dg$lenc / dg$len >= 0.8,]
genes_covered = dg$id
dg = cbind(dg, pi = dg$nd / dg$lenc)
to = ddply(dg, .(fam), summarise, cnt = length(id), q25 = quantile(pi, 0.25), q50 = quantile(pi, 0.5), q75 = quantile(pi, 0.75))
to[order(to$q50, decreasing = T), c(1,2,4)][1:10,]

fo = file.path(diro, "42.genefam.pi.refmap.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

## restrict SNP calling to syntenic regions
diry = file.path(Sys.getenv("misc3"), "comp.vnt")
fa = file.path(diry, "81.cvg.tbl")
ta = read.table(fa, header = F, sep = "\t", as.is = T)
gry = with(ta, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

cat(sum(width(gra)), "\n")
cat(sum(width(gry)), "\n")
cat(sum(width(GenomicRanges::intersect(gra, gry))), "\n")

gra = GenomicRanges::intersect(gra, gry)
snpidxs = intersect_idx1(gsnp, gry)
#gsnp2 = GenomicRanges::intersect(gsnp, gry)

nds = intersect_score(grg, gsnp[snpidxs])
bps = intersect_basepair(grg, gra)

tg2 = cbind(tg, nd = nds, bp = bps)
gb = group_by(tg2, id)
dg = dplyr::summarise(gb, fam = fam[1], len = sum(end-beg+1), lenc = sum(bp), nd = sum(nd))

dg = dg[dg$lenc / dg$len >= 0.8,]
genes_covered = dg$id
dg = cbind(dg, pi = dg$nd / dg$lenc)
to = ddply(dg, .(fam), summarise, cnt = length(id), q25 = quantile(pi, 0.25), q50 = quantile(pi, 0.5), q75 = quantile(pi, 0.75))
to[order(to$q50, decreasing = T), c(1,2,4)][1:10,]

fo = file.path(diro, "42.genefam.pi.refmap.syn.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
