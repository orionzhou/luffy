require(GenomicRanges)
require(plyr)
require(rtracklayer)
require(Rsamtools)
require(ggplot2)
require(grid)
source('Location.R')
source('comp.fun.R')
require(VennDiagram)

orgs = get_orgs()
orgs = get_orgs('ingroup')
orgs = c("HM010", "HM004", "HM023")
chrs = sprintf("chr%s", 1:8)

tt = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
tt = data.frame(chr = tt$V1, beg = 1, end = tt$V2)
tt = tt[tt$chr %in% chrs,]
grt = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

tp = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
colnames(tp) = c("chr", "beg", "end")
tp = tp[tp$chr %in% chrs,]
grp = with(tp, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

grf = GenomicRanges::setdiff(grt, grp)

org = orgs[3]

##### venn diagram
dirm = file.path(Sys.getenv('misc3'), "hapmap/12_ncgr")
fm = sprintf("%s/44_snp/%s.tbl", dirm, org)
tm = read.table(fm, sep = '\t', header = T, as.is = T)
grm1 = with(tm[tm$gt == 1,], GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))
grm2 = with(tm[tm$gt == 2,], GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))

dirv = sprintf("%s/%s_HM101/23_blat", Sys.getenv("misc3"), org)
fv = sprintf("%s/31.9/snp", dirv)
tv = read.table(fv, sep = '\t', header = F, as.is = T)
colnames(tv) = c("chr", "pos", "ref", "alt", "qid", "qpos", "cid", "lev")
tv = tv[tv$lev == 1 & tv$chr %in% chrs, c(1:2,4)]
grv = with(tv, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))


  tc = merge(tm[tm$gt==2,], tv, by = c('chr', 'pos'), all = T)
  tc = cbind(tc, type = NA)
  tc$type[!is.na(tc$alt.x) & !is.na(tc$alt.y)] = 'ovl'
  tc$type[is.na(tc$alt.x) & !is.na(tc$alt.y)] = 'anm'
  tc$type[!is.na(tc$alt.x) & is.na(tc$alt.y)] = 'mna'
  
  t_ovl = tc[!is.na(tc$alt.x) & !is.na(tc$alt.y),]
  t_anm = tc[is.na(tc$alt.x) & !is.na(tc$alt.y),]
  t_mna = tc[!is.na(tc$alt.x) & is.na(tc$alt.y),]
  n_con = sum(t_ovl$alt.x == t_ovl$alt.y)
  n_dis = sum(t_ovl$alt.x != t_ovl$alt.y)

  n_dis / (n_con + n_dis)

gr_ovl = with(t_ovl, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))
gr_mna = with(t_mna, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))
gr_anm = with(t_anm, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))

summary(t_ovl$qual)
summary(t_mna$qual)
summary(t_ovl$rd)
summary(t_mna$rd)

a1 = sum(tm$gt == 2)
a2 = nrow(tv) 
cr = nrow(t_ovl) 
venn.plot <- draw.pairwise.venn(area1 = a1, area2 = a2, cross.area = cr, 
  category = c("Mapping-based Calls", "Synteny-based Calls"),
  fill = c("dodgerblue", "firebrick"), lty = "solid",
  cex = 1, cat.cex = 1, cat.pos = c(200, 20), cat.dist = 0.04,
#  cat.just = list(c(-1, -1), c(1, 1)),
)

# enrichment of mapping-only calls in SV regions
fs = sprintf("%s/%s_HM101/31_sv/11.sv.tbl", Sys.getenv("misc3"), org)
ts = read.table(fs, sep = '\t', header = T, as.is = T)
ts = ts[ts$chr %in% chrs,]
gr_sv = with(ts, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

fx = sprintf("%s/%s_HM101/23_blat/31.9/gax", Sys.getenv("misc3"), org)
tx = read.table(fx, sep = '\t', header = F, as.is = T)
tx = tx[tx$V1 %in% chrs,]
gr_gax = with(tx[tx$V10==1,], GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
gr_var = with(tx[tx$V10>1,], GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))


grv = reduce(c(gr_sv, gr_var))
gru = GenomicRanges::setdiff(grf, reduce(c(gr_gax, grv)))

tpc1 = data.frame(type = 'genome', cla = c('synteny', 'sv', 'uncovered'), 
  cnt = c(sum(width(gr_gax)), sum(width(grv)), sum(width(gru))), 
  total = sum(width(grf)))
tpc2 = data.frame(type = 'mna', cla = c('synteny', 'sv', 'uncovered'), 
  cnt = c(sum(width(GenomicRanges::intersect(gr_mna, gr_gax))), 
  sum(width(GenomicRanges::intersect(gr_mna, grv))), 
  sum(width(GenomicRanges::intersect(gr_mna, gru)))), 
  total = sum(width(gr_mna)))
tpc = rbind(tpc1, tpc2)
tpc = cbind(tpc, pct = tpc$cnt / tpc$total)

p_mna_sv = ggplot(tpc) + 
  geom_bar(aes(x = type, y = pct, fill = cla), position = 'fill', stat = 'identity', width = 0.8) +
#  scale_fill_manual(values = c('mediumseagreen', 'lightsalmon', 'burlywood1')) +
  scale_fill_brewer(palette = "Pastel1", breaks = c("synteny", "sv", "uncovered"), labels = c("regions covered by synteny", "regions affected by SV (resolved)", "other regions")) +
  scale_x_discrete(name = 'proportion', expand = c(0, 0), labels = c("genome-wide total", "Mapping-based SNP calls")) +
  scale_y_continuous(name = '', expand = c(0, 0)) +
  coord_polar(theta = 'y') +
  theme_bw() +
  theme(legend.position = "top", legend.direction = "vertical", legend.title = element_blank(), legend.text = element_text(size = 8), legend.key.size = unit(0.8, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,1,0,0), "lines")) +
  theme(axis.title = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "brown", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 60))
#p_mna_sv

# enrichment of hets in mna sets
  htc = merge(tm[tm$gt<=2,], tv, by = c('chr', 'pos'), all = T)
  
  ht_ovl = htc[!is.na(htc$alt.x) & !is.na(htc$alt.y),]
  ht_anm = htc[is.na(htc$alt.x) & !is.na(htc$alt.y),]
  ht_mna = htc[!is.na(htc$alt.x) & is.na(htc$alt.y),]

tmp1 = table(ht_ovl$gt)
tmp2 = table(ht_mna$gt)
tht = data.frame('snptype' = rep(c('het', 'hom'), 2), 
  'valid' = rep(c('ovl', 'mna'), each = 2),
  cnt = as.numeric(c(tmp1, tmp2)),
  total = rep(c(sum(tm$gt==1), sum(tm$gt==2)), 2))
tht = cbind(tht, pct = tht$cnt / tht$total)

p_het_ovl = ggplot(tht) + 
  geom_bar(aes(x = snptype, y = pct, fill = valid), position = 'fill', stat = 'identity', width = 0.8) +
#  scale_fill_manual(values = c('mediumseagreen', 'lightsalmon', 'burlywood1')) +
  scale_fill_brewer(palette = "Pastel1", breaks = c("ovl", "mna"), labels = c("Validated by synteny-based approach", "Not validated by synteny-based approach")) +
  scale_x_discrete(name = 'proportion', expand = c(0, 0), labels = c("heterozygous calls", "homozygous calls")) +
  scale_y_continuous(name = '', expand = c(0, 0)) +
  coord_polar(theta = 'y') +
  theme_bw() +
  theme(legend.position = "top", legend.direction = "vertical", legend.title = element_blank(), legend.text = element_text(size = 8), legend.key.size = unit(0.8, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(1,1,1,0), "lines")) +
  theme(axis.title = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "brown", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 60))
#p_het_ovl

# enrichment of het calls (mna) in high-divergence regions
x = tt$end
names(x) = tt$chr
gr = tileGenome(x, tilewidth = 1000, cut.last.tile.in.chrom = T)

tw = data.frame(chr = seqnames(gr), beg = start(gr), end = end(gr), 
  len = width(gr), stringsAsFactors = F)

bp_gap = intersect_basepair(gr, grp)
cnt_het = intersect_count(gr, grm1)
cnt_hom = intersect_count(gr, grm2)

to = cbind(tw, len_ng = tw$len - bp_gap, cnt_het = cnt_het, cnt_hom = cnt_hom)
to = to[to$len_ng >= to$len * 0.4,]
to = cbind(to, snpd = (to$cnt_het+to$cnt_hom)/to$len_ng, prop_het = to$cnt_het / (to$cnt_het+to$cnt_hom) * 100)
to = to[!is.na(to$prop_het),]
intvs = c(seq(0, 0.1, 0.005), Inf)
labs = 100 - 100 * intvs[-1]
labs[length(labs)] = "<90.0"
to = cbind(to, snpd2 = cut(to$snpd, intvs))
levs = sort(unique(to$snpd2))

tx = ddply(to, .(snpd2), summarise, q25 = quantile(prop_het, 0.25), q50 = quantile(prop_het, 0.5), q75 = quantile(prop_het, 0.75))

p_het_idty = ggplot(tx) +
  geom_crossbar(aes(x = snpd2, y = q50, ymin = q25, ymax = q75), 
    geom_params = list(width = 0.8)) + 
  scale_x_discrete(name = 'Percent Identify', expand = c(0, 0), breaks = levs[seq(1,41,5)], labels = labs[seq(1,41,5)]) +
  scale_y_continuous(name = '% heterozygous calls', expand = c(0.02, 0)) +
  theme_bw() +
#  theme(legend.position = "top", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,1,0,0), "lines")) +
  theme(axis.title = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 90, hjust = 0.5))
#p_het_idty

# plot rd and qual distribution for mna/ovl calls
intvs = c(seq(0,60,2), Inf)
tmp1 = table(cut(t_ovl$rd, breaks = intvs))
tmp2 = table(cut(t_mna$rd, breaks = intvs))
to1 = data.frame(intv = names(tmp1), cnt = as.numeric(tmp1), type = 'ovl', stringsAsFactors = F)
to1 = cbind(to1, dens = to1$cnt/sum(to1$cnt))
to2 = data.frame(intv = names(tmp2), cnt = as.numeric(tmp2), type = 'mna', stringsAsFactors = F)
to2 = cbind(to2, dens = to2$cnt/sum(to2$cnt))
to = rbind(to1, to2)
to$intv = factor(to$intv, levels = names(tmp1))
levs = sort(unique(to$intv))
labs = intvs[-1]
labs[length(labs)] = "60+"

p_rd = ggplot(to) +
  geom_bar(aes(x = intv, y = dens, fill = type), 
    position = 'dodge', stat = 'identity', geom_params = list(width = 0.9)) + 
  scale_fill_manual(values = c('mediumseagreen', 'lightsalmon'), labels = c("SNPs called only by reference mapping", "SNPs called by both approaches")) +
#  scale_fill_brewer(palette = "Accent") +
  scale_x_discrete(name = 'Read Depth', expand = c(0, 0), breaks = levs[seq(1,41,5)], labels = labs[seq(1,41,5)]) +
  scale_y_continuous(name = 'Density', expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = c(0.5, 0.8), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.6, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_blank(), legend.text = element_text(size = 8, angle = 0)) +
  theme(plot.margin = unit(c(0,1,0,0), "lines")) +
  theme(axis.title = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "brown", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0))
#p_rd


intvs = c(seq(3, 14, 0.5), Inf)
tmp1 = table(cut(log(t_ovl$qual), breaks = intvs))
tmp2 = table(cut(log(t_mna$qual), breaks = intvs))
to1 = data.frame(intv = names(tmp1), cnt = as.numeric(tmp1), type = 'ovl', stringsAsFactors = F)
to1 = cbind(to1, dens = to1$cnt/sum(to1$cnt))
to2 = data.frame(intv = names(tmp2), cnt = as.numeric(tmp2), type = 'mna', stringsAsFactors = F)
to2 = cbind(to2, dens = to2$cnt/sum(to2$cnt))
to = rbind(to1, to2)
to$intv = factor(to$intv, levels = names(tmp1))
p_qual = ggplot(to) +
  geom_bar(aes(x = intv, y = dens, fill = type), 
    position = 'dodge', stat = 'identity', geom_params = list(width = 0.9)) + 
  scale_fill_manual(values = c('mediumseagreen', 'lightsalmon')) +
#  scale_fill_brewer(palette = "Accent") +
  scale_x_discrete(name = 'qual', expand = c(0, 0)) +
  scale_y_continuous(name = '# SNPs', expand = c(0, 0)) +
  theme_bw() +
#  theme(legend.position = "top", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(0,1,1,0), "lines")) +
  theme(axis.title = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "grey", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 90))
#p_qual

# look at sequence percent identity around called SNPs
x = tt$end
names(x) = tt$chr
gr = tileGenome(x, tilewidth = 1000, cut.last.tile.in.chrom = T)

tw = data.frame(chr = seqnames(gr), beg = start(gr), end = end(gr), 
  len = width(gr), stringsAsFactors = F)

bp_gap = intersect_basepair(gr, grp)
bp_gax = intersect_basepair(gr, gr_gax)
cnt_ovl = intersect_count(gr, gr_ovl)
cnt_anm = intersect_count(gr, gr_anm)

to = cbind(tw, len_ng = tw$len - bp_gap, len_gax = bp_gax, cnt_ovl = cnt_ovl, cnt_anm = cnt_anm)
to = to[to$len_gax > to$len * 0.3,]
to = cbind(to, snpd = (to$cnt_ovl+to$cnt_anm)/to$len_gax, pct_anm = to$cnt_anm / (to$cnt_anm+to$cnt_ovl) * 100)
to = to[!is.na(to$pct_anm),]
intvs = c(seq(0, 0.2, 0.005), Inf)
labs = 100 - 100 * intvs[-1]
labs[length(labs)] = "<80.0"
to = cbind(to, snpd2 = cut(to$snpd, intvs))
levs = sort(unique(to$snpd2))

tx = ddply(to, .(snpd2), summarise, q25 = quantile(pct_anm, 0.25), q50 = quantile(pct_anm, 0.5), q75 = quantile(pct_anm, 0.75))

p_pct_idt = ggplot(tx) +
  geom_crossbar(aes(x = snpd2, y = q50, ymin = q25, ymax = q75), 
    geom_params = list(width = 0.7, size = 0.3)) + 
  scale_x_discrete(name = 'Percent Identify', expand = c(0, 0), breaks = levs[seq(1,41,5)], labels = labs[seq(1,41,5)]) +
  scale_y_continuous(name = '% synteny-only calls', expand = c(0.02, 0)) +
  theme_bw() +
#  theme(legend.position = "top", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,1,0,0), "lines")) +
  theme(axis.title = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 90, hjust = 0.5))
#p_pct_idt

### combined plot
fo = sprintf("%s/comp.stat/snp_%s.pdf", Sys.getenv("misc3"), org)
pdf(file = fo, width = 6, height = 8, bg = 'transparent')
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2, heights = c(3,2.5,2.5))))

vlt <- viewport(layout.pos.row = 1, layout.pos.col = 1)
pushViewport(vlt)
grid.draw(venn.plot)
popViewport()

print(p_mna_sv, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_rd, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p_pct_idt, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(p_het_ovl, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(p_het_idty, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))

dco = data.frame(x = rep(1:3, each = 2), y = rep(1:2, 3), lab = LETTERS[1:6])
for (i in 1:nrow(dco)) {
  x = dco$x[i]; y = dco$y[i]; lab = dco$lab[i]
  grid.text(lab, x = 0, y = unit(1, 'npc'), just = c('left', 'top'), gp = gpar(col = "black", fontface = 2, fontsize = 20),
    vp = viewport(layout.pos.row = x, layout.pos.col = y))
}
dev.off()


