require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(gtable)
require(grid)
require(RColorBrewer)
#require(Gviz)
source('Location.R')
source('comp.fun.R')


prepare_data <- function(chr, beg, grt, grp, tg, grl, grc, grn, grvs, grvl) {
beg = beg * 1000000 + 1
winsize = 50000
winnum = 1000000 / winsize
begs = seq(beg, by = winsize, length.out = winnum)
ends = begs + winsize - 1
tw = data.frame(chr = chr, beg = begs, end = ends, len = winsize, stringsAsFactors = F)
gr = with(tw, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
  
bp_gap = intersect_basepair(gr, grp)
bp_nogap = tw$len - bp_gap

tw = cbind(tw, len_ng = bp_nogap)
#tw = tw[tw$chr != 'chrU' & tw$len_ng/tw$len > 0.5,]
gr = with(tw, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

### gene family stats - basepair %
h = list(
  	"TE" = "TE",
  	"NBS_LRR" = c("CC-NBS-LRR", "TIR-NBS-LRR", "NB-ARC", "TIR"),
  	"RLK" = c("LRR-RLK", "RLK"),
  	'NCR' = 'NCR',
  	'LRR' = 'LRR',
  	"F_box" = "F-box"
)

dg = tw[,1:3]
for (key in names(h)) {
  tgs = tg[tg$cat %in% h[[key]],]
  gs = with(tgs, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
  gd = intersect_basepair(gr, reduce(gs)) / tw$len_ng
  dg = cbind(dg, gd = gd)
  colnames(dg)[length(dg)] = key
}

### proportion covered based in synteny
dy = tw[,1:3]
for (qname in qnames_15) {
  len_syn = intersect_basepair(gr, grl[[qname]])
  pc = len_syn / tw$len_ng
  pc[tw$len_ng == 0] = 0
  dy = cbind(dy, org = pc)
  colnames(dy)[length(dy)] = qname
}

### theta-pi, covered bases in 9/12 accessions
lenc = intersect_basepair(gr, grc)
pcov = lenc / tw$len_ng

nds = intersect_score(gr, grn)
pi_snp = nds / lenc
ndi = intersect_score(gr, grvs)
pi_indel = ndi / lenc
ndv = intersect_score(gr, grvl)
pi_sv = ndv / lenc

pi_snp[pcov < 0.2] = NA
pi_indel[pcov < 0.2] = NA
pi_sv[pcov < 0.2] = NA
ds = cbind(tw[,1:3], Gaps = 1-tw$len_ng/tw$len, Coverage = lenc/tw$len_ng, Pi_SNP = pi_snp, Pi_InDel = pi_indel, Pi_SV = pi_sv)

  list(tw = tw, dg = dg, dy = dy, ds = ds)
}

sub_plots <- function(chr, beg, tw, dg, dy, ds) {

ptitle = sprintf("%s:%02dMb", chr, beg)
##### sliding window plot (100 x 10kb)
dyw = cbind(dy[,c(-1:-3)], beg = dy[,2])
dyl = reshape(dyw, direction = 'long', varying = list(1:15), idvar = c('beg'), timevar = 'org', v.names = 'syn', times = colnames(dyw)[1:15])
dyl$org = factor(dyl$org, levels = rev(qnames_15))

cols1 = rev(c("#282b68", "#324387", "#385193", "#498bbd", "#71c5cd", "#81c185", "#afcf45", "#faed29", "#ea862d", "#db382b", "#bb242a"))
colsy = brewer.pal(11, "Spectral")

vmin = min(dyl$syn, na.rm = T)
vmax = max(dyl$syn, na.rm = T)
breaks = seq(vmin, vmax, length.out = length(colsy) + 1)
dyl = cbind(dyl, syn2 = cut(dyl$syn, breaks, include.lowest = T))

labs = rep('', length(colsy))
ltitle = sprintf("%.03f-%.03f", vmin, vmax)

p_syn <- ggplot(dyl) +
  geom_tile(aes(x = beg, y = org, fill = syn2, height = 0.8)) +
  theme_bw() + 
  scale_x_continuous(name = '', expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0), name = '') +
  scale_fill_manual(name = ltitle, labels = labs, values = colsy, guide = guide_legend(nrow = 1, byrow = T, label = T, label.position = 'top')) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.title = element_text(size = 7), legend.key.size = unit(0.3, 'lines'), legend.key.width = unit(0.3, 'lines'), legend.text = element_blank(), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "cm")) +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, color='black', size=1)) + #deeppink4
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0,'cm')) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))

#fo = file.path(dirw, "33.win.1.pdf")
#ggsave(py, filename = fo, width = 4, height = 4)


## different statistics (gene density + gap + coverage + theta-pi)
dw = cbind(ds, dg[,-1:-3])
opts = colnames(dw)[-1:-3]
dcl = data.frame()
for (opt in opts) {
	cols = c("white", brewer.pal(9, "Reds"))
	vals = dw[,opt]
	vmin = min(dw[,opt], na.rm = T); vmax = max(dw[,opt], na.rm = T)
	breaks = seq(vmin, vmax, length.out = length(cols)+1)
	if(vmin == vmax) {
		itvs = rep(vmin, length(vals))
		itvs = factor(itvs, levels = seq(vmin, by = 0.1, length.out = length(cols)))
	} else {
		itvs = cut(vals, breaks, include.lowest = T)
	}

	#levs = sort(unique(itvs))
	#nlev = length(levs)
	#itvs = factor(itvs, levels = levs)
	#if(nlev == 1) {
	#	cols = 'white'
	#} else if(nlev == 2) {
	#	cols = c('white', 'black')
	#} else if(nlev < 9) {
	#	cols = brewer.pal(max(3,nlev), "Greys")
	#}
	col_map = cols
	names(col_map) = levels(itvs)
	vals_mapped = as.character(col_map[itvs])
	vals_mapped[is.na(vals_mapped)] = 'grey75'
	
	labg = rep('', length(cols))
	ltitle = sprintf("%% bases (%.01f-%.01f)", vmin, vmax)
	dcl1 = cbind(tw[,1:2], opt = opt, col = vals_mapped)
	dcl = rbind(dcl, dcl1)
}
dcl$opt = factor(dcl$opt, levels = rev(opts))

p_sta <- ggplot(dcl) +
  geom_tile(aes(x = beg, y = opt, fill = col, height = 0.8)) +
  theme_bw() + 
  scale_x_continuous(name = ptitle, expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0), name = '') +
  scale_fill_identity() +
  theme(legend.position = 'none') + 
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, color='black', size=1)) + #deeppink4
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_text(colour = "blue", size = 8), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))

	list(syn = p_syn, sta = p_sta)
}
