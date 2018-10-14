require(rtracklayer)
require(xlsx)
source("comp.fun.R")
source("comp.plot.fun.R")

dirw = file.path(Sys.getenv("misc3"), 'comp.stat', 'figs')
fl = file.path(dirw, 'loci.xlsx')

##### fine-scale synteny plot
source("comp.plot.fun.R")
tl = read.xlsx(fl, sheetIndex = 1, header = T)

tracks = c('tgene', 'taxis', 'tgap', 'link', 'qgap', 'qaxis', 'qgene')
#tracks = c('taxis', 'tgap', 'link', 'qgap', 'qaxis')
i = 28
tls = tl[tl$i == i,]

gro =  with(tls, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
qnames = strsplit(as.character(tls$qnames), split = ' ')[[1]]

cfgs = get_genome_cfgs(c(tname, qnames))
dats = prep_plot_data(gro, cfgs, tname, qnames, tracks)
res = comp.plot(dats, tname, qnames, tracks, draw.legend.gene = T, scale.ht  = unit(0.8, 'npc'))

fn = sprintf("%s/fig%03d.pdf", dirw, i)
CairoPDF(file = fn, width = 7, height = res$ht/72, bg = 'transparent')
grid.newpage()
grid.draw(res$grobs)
dev.off()

#### large-scale synteny plot
source("comp.plot.fun.R")
tl = read.xlsx(fl, sheetIndex = 1, header = T)

tracks = c('taxis', 'link', 'qaxis')
tracks = c('taxis', 'tgap', 'link', 'qgap', 'qaxis')
i = 93
tls = tl[tl$i == i,]

gro =  with(tls, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
qnames = strsplit(as.character(tls$qnames), split = ' ')[[1]]

cfgs = get_genome_cfgs(c(tname, qnames))
dats = prep_plot_data(gro, cfgs, tname, qnames, tracks, largescale = T)
res = comp.plot(dats, tname, qnames, tracks, scale.ht  = unit(0.85, 'npc'), largescale = T)

fn = sprintf("%s/fig%03d.pdf", dirw, i)
CairoPDF(file = fn, width = 8, height = res$ht/72, bg = 'transparent')
grid.newpage()
grid.draw(res$grobs)
dev.off()

#### chr4-8 translocation
source("comp.plot.fun.R")
tl = read.xlsx(fl, sheetIndex = 1, header = T)

tracks = c('taxis', 'link', 'qaxis')
i = 92
tls = tl[tl$i == i,]

gro =  with(tls, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
qnames = strsplit(as.character(tls$qnames), split = ' ')[[1]]

cfgs = get_genome_cfgs(c(tname, qnames))
qods = list("HM340.FN" = c("scf010", "scf005", "scf015", "scf002"))
dats = prep_plot_data(gro, cfgs, tname, qnames, tracks, largescale = T, fill_prop = 0.05, qods = qods)
res = comp.plot(dats, tname, qnames, tracks, scale.ht  = unit(0.45, 'npc'), largescale = T)

fn = sprintf("%s/fig%03d.pdf", dirw, i)
CairoPDF(file = fn, width = 7, height = res$ht/72, bg = 'transparent')
grid.newpage()
grid.draw(res$grobs)

#xsu = c(165, 328, 333, 340)
#ysu = c(1:3) * 105 - 30
#xs = rep(xsu, length(ysu))
#ys = rep(ysu, each = length(xsu))
#grid.text(rep("x", length(xs)), x = unit(xs, 'points'), y = unit(ys, 'points'), 
#  just = c('left', 'bottom'), 
#  gp = gpar(col = "navy", fontface = 1, fontsize = 8)
#)
dev.off()

## zoom in
source("comp.plot.fun.R")
tl = read.xlsx(fl, sheetIndex = 1, header = T)

tracks = c('taxis', 'link', 'qaxis')
i = 93
tls = tl[tl$i == i,]

gro =  with(tls, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
qnames = strsplit(as.character(tls$qnames), split = ' ')[[1]]

cfgs = get_genome_cfgs(c(tname, qnames))
qods = list("HM340.FN" = c("scf010", "scf005", "scf015", "scf002"))
dats = prep_plot_data(gro, cfgs, tname, qnames, tracks, largescale = T, fill_prop = 0.01, qods = qods)
res = comp.plot(dats, tname, qnames, tracks, scale.ht  = unit(0.15, 'npc'), largescale = T)

fn = sprintf("%s/fig%03d.pdf", dirw, i)
CairoPDF(file = fn, width = 7, height = res$ht/72, bg = 'transparent')
grid.newpage()
grid.draw(res$grobs)
dev.off()

## zoom in
source("comp.plot.fun.R")
tl = read.xlsx(fl, sheetIndex = 1, header = T)

tracks = c('taxis', 'link', 'qaxis')
i = 94
tls = tl[tl$i == i,]

gro =  with(tls, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
qnames = strsplit(as.character(tls$qnames), split = ' ')[[1]]

cfgs = get_genome_cfgs(c(tname, qnames))
dats = prep_plot_data(gro, cfgs, tname, qnames, tracks, largescale = F)
res = comp.plot(dats, tname, qnames, tracks, scale.ht  = unit(0.1, 'npc'), largescale = F)

fn = sprintf("%s/fig%03d.pdf", dirw, i)
CairoPDF(file = fn, width = 7, height = res$ht/72, bg = 'transparent')
grid.newpage()
grid.draw(res$grobs)

xsu = c(141, 286, 331, 409)
ysu = c(1:3) * 135 - 45
xs = rep(xsu, length(ysu))
ys = rep(ysu, each = length(xsu))
grid.text(rep("x", length(xs)), x = unit(xs, 'points'), y = unit(ys, 'points'), 
  just = c('left', 'bottom'), 
  gp = gpar(col = "navy", fontface = 1, fontsize = 8)
)
dev.off()

### tandem duplication illustration
chr = "chr8"
beg = 3370714
end = 3416737
gro =  GRanges(seqnames = chr, ranges = IRanges(beg, end = end))
source("comp.plot.fun.R")

qnames = sprintf("HM%03d", c(56, 10, 129, 34))
qnames = c("HM034")
tracks = c('tgene', 'trnaseq', 'taxis', 'tgap', 'link', 'qgap', 'qaxis', 'qgene', 'qpacbio', 'qrnaseq')
  
dats = prep_plot_data(gro, cfgs, tname, qnames, tracks)
res = comp.plot(dats, tname, qnames, tracks, scale.ht = unit(0.8, 'npc'), draw.legend.gene = T)

fn = sprintf("%s/illus_tandup.pdf", dirw)
CairoPDF(file = fn, width = 7, height = res$ht/72, bg = 'transparent')
grid.newpage()
grid.draw(res$grobs)
dev.off()

### SNP calling inconsistency illustration
source("comp.plot.fun.R")
chr = "chr5"
beg = 1555000
end = 1605000
gro =  GRanges(seqnames = chr, ranges = IRanges(beg, end = end))

qnames = c("HM125")
tracks = c('qgene', 'qaxis', 'qgap', 'link', 'tgap', 'taxis', 'tgene', 'tmapp', 'mcov', 'msnp', 'tsnp')

dats = prep_plot_data(gro, cfgs, tname, qnames, tracks)
res1 = comp.plot(dats, tname, qnames, tracks, scale.ht = unit(0.9, 'npc'))

chr = "chr3"
beg = 2900000
end = 2925000
gro =  GRanges(seqnames = chr, ranges = IRanges(beg, end = end))
source("comp.plot.fun.R")

qnames = c("HM034")
tracks = c('qgene', 'qaxis', 'qgap', 'link', 'tgap', 'taxis', 'tgene', 'tmapp', 'mcov', 'msnp', 'tsnp', 'tpct')
  
dats = prep_plot_data(gro, cfgs, tname, qnames, tracks)
res2 = comp.plot(dats, tname, qnames, tracks, scale.ht = unit(0.9, 'npc'))

numrow = 4; numcol = 3
ress = list(res1, res2)
htt = 15; htm = 5; wdm = 5
ht = res1$ht + res2$ht + htt + htm

fn = sprintf("%s/illus_snpcall.pdf", dirw)
CairoPDF(file = fn, width = 7.15, height = ht/72, bg = 'transparent')
grid.newpage()

pushViewport(viewport(layout = grid.layout(numrow, numcol, 
  heights = c(htt, res1$ht, htm, res2$ht),
  widths = unit.c(unit(wdm, 'points'), unit(1, 'npc') - unit(wdm*2, 'points'), unit(wdm, 'points')))
))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(plot_legend_gene())
popViewport()
for (i in 1:length(ress)) {
  pushViewport(viewport(layout.pos.row = 1+i*2-1, layout.pos.col = 2))
  grid.draw(ress[[i]]$grobs)
  
  if(i == 1) {
xs = c(88, 250, 407)
wds = c(155, 63, 35)
ys = rep(3, 3)
hts = rep(103, 3)
grid.rect(x = unit(xs, 'points'), width = unit(wds, 'points'), 
  y = unit(ys, 'points'), height = unit(hts, 'points'),
  just = c('left', 'bottom'),
  gp = gpar(lwd = 1, col = 'dodgerblue3', fill = NA)
)
labs = LETTERS[1:length(xs)]
grid.text(labs, x = unit(xs, 'points'), y = unit(ys+hts, 'points'), 
  just = c('left', 'top'), 
  gp = gpar(col = "dodgerblue3", fontface = 2, fontsize = 13)
)
  } else {
xs = c(390)
wds = c(80)
ys = rep(1, length(xs))
hts = rep(123, length(xs))
grid.rect(x = unit(xs, 'points'), width = unit(wds, 'points'), 
  y = unit(ys, 'points'), height = unit(hts, 'points'),
  just = c('left', 'bottom'),
  gp = gpar(lwd = 1, col = 'dodgerblue3', fill = NA)
)
labs = LETTERS[4:(3+length(xs))]
grid.text(labs, x = unit(xs, 'points'), y = unit(ys+hts, 'points'), 
  just = c('left', 'top'), 
  gp = gpar(col = "dodgerblue3", fontface = 2, fontsize = 13)
)
  }
  popViewport()
}
dev.off()

##### SV illustration
### look for del/ins/inv/tlc
fb = '/home/youngn/zhoup/Data/misc3/HM034_HM101/31_sv/05.stb'
tb = read.table(fb, header = T, sep = "\t", as.is = T)
fx = '/home/youngn/zhoup/Data/misc3/HM034_HM101/31_sv/05.stx'
tx = read.table(fx, header = T, sep = "\t", as.is = T)

# ins/del
td = tx[tx$type == "DEL",]
td[td$tend - td$tbeg > 5000,][21:40,]
tn = tx[tx$type == "INS",]
tn[tn$qend - tn$qbeg > 3000,][1:20,]

# tlc
tc = tx[tx$type %in% c('TLC:INS', 'TLC:DEL'),]
tb = tb[,c(1:4,8:10)]
to = merge(tc, tb, by = 'id')
to = to[order(to$tchr.x, to$tbeg.x),]
to[to$tend.x-to$tbeg.x>10000,][41:60,]

# look for inv
idxs = c()
for (i in 1:(nrow(to)-1)) {
  if(to$id[i] == to$id[i+1] & to$tbeg.x[i] == to$tbeg.x[i+1]) {
    idxs = c(idxs, i)
  }
}


## simple SV
tl = read.xlsx(fl, sheetIndex = 2, header = T)
fn = sprintf("%s/illus_sv.pdf", dirw)
scale.hts = c(0.8, 0.2, 0.26, 0.7)

## NBS SV
tl = read.xlsx(fl, sheetIndex = 3, header = T)
fn = sprintf("%s/illus_sv_nbs.pdf", dirw)
scale.hts = c(0.9, 0.8, 0.85)

## CRP SV
tl = read.xlsx(fl, sheetIndex = 4, header = T)
fn = sprintf("%s/illus_sv_crp.pdf", dirw)
scale.hts = c(0.9, 0.8, 0.85, 0.2, 0.2)

## ALPACA loci selection
tl = read.xlsx(fl, sheetIndex = 5, header = T)
fn = sprintf("%s/illus_alpaca.pdf", dirw)
scale.hts = c(0.25, 0.25, 0.25)

## plotting
source("comp.plot.fun.R")
idxs = sort(unique(tl$i))
tracks = c('tgene', 'taxis', 'tgap', 'link', 'qgap', 'qaxis', 'qgene')
#tracks = c('taxis', 'link','qaxis')
ress = list()
hts = c()
for (i in 1:length(idxs)) {
  idx = idxs[i]
  tls = tl[tl$i == idx,]
  gro =  GRanges(seqnames = tls$chr, ranges = IRanges(tls$beg, end = tls$end))
  qnames = strsplit(as.character(tls$qnames), split = ' ')[[1]]
  cfgs = get_genome_cfgs(c(tname, qnames))


  dats = prep_plot_data(gro, cfgs, tname, qnames, tracks)
  res = comp.plot(dats, tname, qnames, tracks, scale.ht = unit(scale.hts[i], 'npc'))
  ress[[i]] = res
  hts = c(hts, res$ht)
}

nplot = length(idxs)
numrow = nplot*2+1; numcol = 3

htt = 15; htm = 5; wdm = 15

CairoPDF(file = fn, width = 7, height = (htt+sum(hts)+htm*nplot)/72, bg = 'transparent')
grid.newpage()

pushViewport(viewport(layout = grid.layout(numrow, numcol, 
  heights = c(htt, as.vector(rbind(rep(htm, nplot), hts))),
  widths = unit.c(unit(wdm, 'points'), unit(1, 'npc') - unit(wdm*2, 'points'), unit(wdm, 'points')))
))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(plot_legend_gene())
popViewport()
for (i in 1:length(idxs)) {
  pushViewport(viewport(layout.pos.row = 1+i*2, layout.pos.col = 2))
  grid.draw(ress[[i]]$grobs)
  popViewport()
}
labs = LETTERS[1:nplot]
grid.text(labs, x = unit(2, 'points'), 
  y = unit(c(nplot:1) * htm + rev(cumsum(rev(hts))), 'points'), 
  just = c('left', 'top'), 
  gp = gpar(col = "black", fontface = 2, fontsize = 16)
)
dev.off()

#### NBS SV illustration 1
qnames = c("HM058", 'HM129', 'HM004')
gro = GRanges(seqnames = c('chr2'), ranges = IRanges(25741000, end = 25746500))

tracks = c('tgene', 'taxis', 'tgap', 'link', 'qgap', 'qaxis', 'qgene')

dats = prep_plot_data(gro, cfgs, tname, qnames, tracks)
res = comp.plot(dats, tname, qnames, tracks, draw.legend.gene = T, scale.ht  = unit(0.85, 'npc'))

fn = sprintf("%s/illus_nbs.pdf", dirw)
CairoPDF(file = fn, width = 7, height = res$ht/72, bg = 'transparent')
grid.newpage()
grid.draw(res$grobs)
dev.off()

