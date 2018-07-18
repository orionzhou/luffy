require(plyr)
require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(RColorBrewer)
require(gtable)
require(dplyr)
source("comp.fun.R")

dirw = file.path(Sys.getenv('misc3'), 'comp.panseq')
diro = file.path(Sys.getenv("misc3"), "comp.stat")


tlen = read.table(cfgs[[tname]]$size, sep = "\t", as.is = T, header = F)
tgap = read.table(cfgs[[tname]]$gap, sep = "\t", as.is = T, header = F)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, V2)))
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, V3)))
grs = setdiff(grt, grp)

for (org in qnames_12) {
  fgax = sprintf("%s/%s_HM101/23_blat/31.9/gax", Sys.getenv('misc3'), org)
  tgax = read.table(fgax, sep = "\t", header = F, as.is = T)[,1:3]
  gr = with(tgax, GRanges(seqnames = V1, ranges = IRanges(V2, V3)))
  grs = c(grs, gr)
}

gro = as(coverage(grs), "GRanges")
scores = mcols(gro)$score
lens = width(gro)
x = data.frame(n_org = scores, len = lens, stringsAsFactors = F)
tx = ddply(x, .(n_org), summarise, size = sum(len), org = 'mixed')
tx = tx[tx$n_org > 0,]
tx$org[tx$n_org == 1] = "HM101"


##### pan-genome AFS
fi = file.path(dirw, '31.refined.tbl')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti = ti[ti$org %in% qnames_12,]
#ddply(ti, .(org), summarise, len = sum(end - beg + 1))

gb = group_by(ti, cid)
tu1 = dplyr::summarise(gb, n_org = n(), 
  orgs = paste(sort(unique(as.character(org))), collapse = "_"),
  size = sum(end - beg + 1))

orgs = c(tname, qnames_12)
tu2 = ddply(tu1, .(n_org), summarise, size = sum(size))
dt1 = data.frame(n_org = tu2$n_org, size = as.integer(tu2$size/tu2$n_org), org = 'mixed', stringsAsFactors = F)
dt1 = dt1[dt1$n_org != 1,]
dt1 = rbind(dt1, tx[tx$n_org != 1,])
dt1 = ddply(dt1, .(n_org), summarise, size = sum(size), org = org[1])

tus = tu1[tu1$n_org == 1,]
dtt = ddply(tus, .(orgs), summarise, size = sum(size))
dt2 = data.frame(n_org = 1, size = dtt$size, org = dtt$orgs, stringsAsFactors = F)
dt2 = rbind(dt2, tx[tx$n_org == 1,])
dt2$org = factor(dt2$org, levels = orgs)
dt2 = dt2[order(dt2$org),]

to = rbind(dt2, dt1)
to$size = to$size / 1000000

cols = c(brewer.pal(12, 'Set3'), brewer.pal(3, 'Set1')[1], 'gray30')
labs = orgs

to$org = factor(to$org, levels = c(orgs, 'mixed'))
to$n_org = factor(to$n_org, levels = sort(as.numeric(unique(to$n_org))))
to = to[order(to$org, decreasing = T),]
p1 = ggplot(to, aes(x = n_org, y = size, fill = org)) +
  geom_bar(stat = 'identity', position = "stack", width=0.5) +
  scale_fill_manual(name = "Accession-Specific", breaks = labs, labels = labs, values = cols) +
  scale_x_discrete(name = '# Sharing Accessions') +
  scale_y_continuous(name = 'Sequences (Mbp)', expand = c(0, 0), limits = c(0, 260)) +
  theme_bw() +
  ggtitle("A") +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(legend.position = c(0.25, 0.8), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.5, 'lines'), legend.title = element_text(size = 8, angle = 0), legend.text = element_text(size = 8, angle = 0)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))

fp = sprintf("%s/61.pan.genome.afs.pdf", diro)
ggsave(p1, filename = fp, width = 4, height = 5)

##### pan/core-genome size curve
fi = file.path(dirw, '31.refined.tbl')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti = within(ti, {len = end - beg + 1})
ti = ti[ti$len > 20,] ##### filter very short segments

tcfg = get_genome_cfg(tname)
tt1 = read.table(tcfg$size, sep = "\t", as.is = T)
grt1 = with(tt1, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt2 = read.table(tcfg$gap, sep = "\t", as.is = T)
grt2 = with(tt2, GRanges(seqnames = V1, ranges = IRanges(V2 + 1, end = V3)))
grt = setdiff(grt1, grt2)

grgs = list()
qnames = qnames_12
for (qname in qnames) {
  fgax = sprintf("%s/%s_%s/23_blat/31.9/gax", Sys.getenv("misc3"), qname, tname)
  gax = read.table(fgax, header = F, sep = "\t", as.is = T)
  grg = with(gax, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
  grgs[[qname]] = grg
}

reps = 1:20
n_orgs = 1:(1+length(qnames))
tp = data.frame(rep = rep(reps, each = length(n_orgs)), 
  n_org = rep(n_orgs, length(rep)), core = NA, pan = NA)
  
for (rep in reps) {
  grc = grt
  grp = grt
  core = sum(width(grc))
  pan = sum(width(grc))
  tp$core[tp$rep == rep && tp$n_org == 1] = core
  tp$pan[tp$rep == rep && tp$n_org == 1] = pan
  
  set.seed(rep*100)
  qnames.rep = sample(qnames)
  for (i in 1:length(qnames.rep)) {
    qname = qnames.rep[i]
    org_str = paste(c(tname, qnames.rep[1:i]), collapse = "+")
    
    grc = intersect(grc, grgs[[qname]])
    core = sum(width(grc))
    tp$core[tp$rep == rep & tp$n_org == i+1] = core
    
    tis = ti[ti$org %in% qnames.rep[1:i], ]
#    tu1 = ddply(tis, .(cid), summarise, n_org = length(org), size = sum(len), .parallel = T)
    tis_df = group_by(tis, cid)
    tu1 = summarise(tis_df, n_org = length(org), size = sum(len))
#    system.time(summarise(tis_df, n_org = length(org), size = sum(len)))
    tu2 = ddply(tu1, .(n_org), summarise, size = sum(size))
    tu3 = cbind(tu2, persize = tu2$size / tu2$n_org)
    pan = sum(width(grp)) + sum(tu3$persize)
    tp$pan[tp$rep == rep & tp$n_org == i+1] = pan
    
    cat(rep, qname, core, pan, "\n")
  }
}
fo = file.path(diro, "63.pan.genome.size.tbl")
write.table(tp, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

##### combo-plot
fi = file.path(diro, "63.pan.genome.size.tbl")
tp = read.table(fi, header = T, sep = "\t", as.is = T)
tp$rep = factor(tp$rep, levels = 1:max(tp$rep))
tp$pan = tp$pan / 1000000
tp$core = tp$core / 1000000

## loess fitting + plot
tp2 = ddply(tp, .(n_org), summarise, panl = quantile(pan, 0.05), panu = quantile(pan, 0.95), corel = quantile(core, 0.05), coreu = quantile(core, 0.95))
p2 = ggplot(tp) +
#  geom_point(aes(x = n_org, y = pan), shape = 1, size = 0.3) +
#  geom_point(aes(x = n_org, y = core), shape = 4, size = 0.3) +
#  geom_boxplot(aes(x = n_org, y = pan, group = n_org), width = 0.4, lwd = 0.2, outlier.color = NA, outlier.size = 0) +
#  geom_boxplot(aes(x = n_org, y = core, group = n_org), width = 0.4, lwd = 0.2, outlier.color = NA, outlier.size = 0) +
  geom_errorbar(data = tp2, mapping = aes(x = n_org, ymin = panl, ymax = panu), stat = 'identity', width = 0.3, lwd = 0.3) +
  geom_errorbar(data = tp2, mapping = aes(x = n_org, ymin = corel, ymax = coreu), stat = 'identity', width = 0.3, lwd = 0.3) +
#  geom_text(aes(x = n_org, y = 0, label = org), geom_params=list(size = 2.5, vjust = 0, angle = 30)) +
  stat_smooth(aes(x = n_org, y = pan, col = 'a'), fill = 'azure4', size = 0.3, se = F) +
  stat_smooth(aes(x = n_org, y = core, col = 'b'), fill = 'azure4', size = 0.3, se = F) +
#  scale_shape(name = "", solid = FALSE, guide = F) +
  scale_color_manual(name = "", labels = c('Pan-genome', 'Core-genome'), values = c("firebrick1", "dodgerblue")) +
  scale_x_continuous(name = '# Genomes Sequenced') +
  scale_y_continuous(name = 'Genome size (Mbp)', expand = c(0, 0), limits = c(0, 480)) + 
  theme_bw() +
  ggtitle("B") +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(legend.position = c(0.25, 0.2), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA), legend.key.size = unit(1, 'lines'), legend.title = element_blank(), legend.text = element_text(size = 8, angle = 0)) +
#    theme(panel.grid = element_blank()) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, color = "blue")) +
  theme(axis.text.y = element_text(size = 8, color = "brown", angle = 90, hjust  = 0.5))

fp = sprintf("%s/63.pan.genome.size.pdf", diro)
ggsave(p2, filename = fp, width = 4, height = 5)


## exponential fitting + extrapolating plot
require(poweRlaw)
tt = ddply(tp, .(n_org), summarise, 
  pan0 = min(pan), pan25 = quantile(pan, 0.25), pan50 = median(pan), pan75 = quantile(pan, 0.75), pan100 = max(pan), panm = mean(pan), 
  core0 = min(core), core25 = quantile(core, 0.25), core50 = median(core), core75 = quantile(core, 0.75), core100 = max(core), corem = mean(core))

#z <- glm(panm ~ b0 + b1 * n_org^b2, data = tt)
# y = Ax^B + C   Ae^Bx + C  b0 + b1*(1-exp(-exp(lrc) * x))
xs = 1:13
xsp = 1:50

#fitl <- loess(pan ~ n_org, tp)
#ysl = predict(fitl, data.frame(n_org = xs), se = F)
dft = data.frame(x = tt$n_org, y = tt$panm)
fite = NLSstAsymptotic(sortedXyData(expression(x),expression(y), dft))
b0 = fite['b0']; b1 = fite['b1']; lrc = fite['lrc']
pansf = b0 + b1*(1-exp(-exp(lrc) * xsp))
b0+b1

#fitl <- loess(core ~ n_org, tp)
#ysl = predict(fitl, data.frame(n_org = xs), se = F)
dft = data.frame(x = tt$n_org, y = tt$corem)
fite = NLSstAsymptotic(sortedXyData(expression(x),expression(y), dft))
b0 = fite['b0']; b1 = fite['b1']; lrc = fite['lrc']
coresf = b0 + b1*(1-exp(-exp(lrc) * xsp))
b0+b1

tf = data.frame(x=xsp, panf = pansf, coref = coresf)

tt$n_org = factor(tt$n_org)
p2 = ggplot(tp) +
  geom_linerange(data = tt, mapping = aes(x = n_org, y = pan50, ymin = pan0, ymax = pan100), col = 'darkorchid4', size = 0.5) +
  geom_linerange(data = tt, mapping = aes(x = n_org, y = core50, ymin = core0, ymax = core100), col = 'darkorchid4', size = 0.5) +
  geom_jitter(aes(x = n_org, y = pan), size = 0.2, width = 0.5) +
  geom_jitter(aes(x = n_org, y = core), size = 0.2, width = 0.5) +
  geom_line(data = tf, aes(x = x, y = panf, col = 'a'), size = 0.5) +
  geom_line(data = tf, aes(x = x, y = coref, col = 'b'), size = 0.5) +
  scale_color_manual(name = "", labels = c('Pan-genome', 'Core-genome'), values = c("firebrick1", "dodgerblue")) +
  scale_x_discrete(name = '# Genomes Sequenced') +
  scale_y_continuous(name = 'Genome size (Mbp)', expand = c(0, 0), limits = c(220, 450)) + 
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(legend.position = c(0.85, 0.65), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA), legend.key.size = unit(1, 'lines'), legend.title = element_blank(), legend.text = element_text(size = 8, angle = 0)) +
  theme(plot.margin = unit(c(1,1,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, color = "blue")) +
  theme(axis.text.y = element_text(size = 8, color = "grey", angle = 90, hjust  = 0.5))

fp = sprintf("%s/63.pan.genome.size.pdf", diro)
ggsave(p2, filename = fp, width = 5, height = 4)

### 2-column plot
g1 = ggplotGrob(p1)
g2 = ggplotGrob(p2)
gs = list(g1, g2)
wids = c(4,4)
g <- gtable_matrix(name = 'demo', grobs = matrix(gs, ncol = length(gs)), widths = wids, heights = 1)

fo = sprintf("%s/63.pan.genome.pdf", diro)
pdf(file = fo, width = 8, height = 5, bg = 'transparent')
grid.newpage()
grid.draw(g)
dev.off()


