require(plyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(RColorBrewer)
require(gtable)
source('comp.fun.R')

diro = file.path(Sys.getenv("misc3"), "comp.stat")

##### compare InDels
chrs = sprintf("chr%s", 1:8)

dirm = file.path(Sys.getenv('misc3'), "hapmap/12_ncgr")
fm = sprintf("%s/52.stat.tbl", dirm)
tm = read.table(fm, sep = '\t', header = T, as.is = T)
tms = tm[tm$chr %in% chrs & (tm$rsize == 1 | tm$asize == 1) & tm$rsize <= 50 & tm$asize <= 50,]

dirv = sprintf("%s/comp.vnt", Sys.getenv("misc3"))
fv = sprintf("%s/52.stat.tbl", dirv)
tv = read.table(fv, sep = '\t', header = T, as.is = T)
tv = cbind(tv, size = (tv$rsize + tv$asize - 2))
tvs = tv[tv$chr %in% chrs & (tv$rsize == 1 | tv$asize == 1) & tv$rsize <= 50 & tv$asize <= 50 & tv$rsize != tv$asize,]

x = c(-49:-1, 1:49)
  
tb = table(tms$asize - tms$rsize)
y = tb[as.character(x)]
y[is.na(y)] = 0
df1 = data.frame(type = 'Mapping-based', len = x, cnt = y)

tb = table(tvs$asize - tvs$rsize)
y = tb[as.character(x)]
y[is.na(y)] = 0
df2 = data.frame(type = 'Assembly-based', len = x, cnt = y)

do = cbind(rbind(df1, df2))
do$cnt = log10(do$cnt)
do$cnt[!is.finite(do$cnt)] = 0

p = ggplot(do) +
  geom_bar(mapping = aes(x = len, y = cnt, fill = type), 
    stat = 'identity', position = 'dodge', 
    geom_params=list(width = 0.8, alpha = 0.8)) +
  scale_fill_brewer(palette='Set1', name = '', labels = c("Reference mapping-based calls", "Synteny-base calls")) +
  scale_x_continuous(name = 'Indel Size (bp)', expand = c(0, 0), breaks = seq(-40, 40, by = 10)) +
  scale_y_continuous(name = '# events (log10)', limits = c(0, 6), expand = c(0, 0), breaks = 0:6)+#, labels = expression(10[])) +
  theme_bw() +
  theme(legend.position = "top", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(0,1,1,0), "lines")) +
  theme(plot.margin = unit(c(0.5,1,0,0), "lines")) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0))

fo = file.path(diro, "comp_idm_size.pdf")
ggsave(p, filename = fo, width = 5, height = 4)


do = data.frame()
for (org in c("HM058", "HM060", "HM095", "HM034")) {
#org = "HM034"
  dira = sprintf("%s/%s_HM101/23_blat/31.9", Sys.getenv('misc3'), org)
  fa = file.path(dira, "vnt.tbl")
  ta = read.table(fa, sep = "\t", header = F, as.is = T)
  colnames(ta) = c("chr", "pos", "ref", "alt", "score")

  dirm = file.path(Sys.getenv("misc3"), "hapmap", "12_ncgr", "44_tbl")
  fm = sprintf("%s/%s.tbl", dirm, org)
  tm = read.table(fm, sep = "\t", header = F, as.is = T)
  colnames(tm) = c("chr", "pos", "ref", "alt", "score")

  tms = tm[tm$chr %in% chrs & (nchar(tm$ref) != 1 | nchar(tm$alt) != 1), ]
  tas = ta[ta$chr %in% chrs & (nchar(ta$ref) != 1 | nchar(ta$alt) != 1), ]
#  tc = merge(tas[,1:4], tms, by = c('chr', 'pos'), all = T)
  
  x = c(-49:-1, 1:49)
  
  tb = table(nchar(tms$alt) - nchar(tms$ref))
  y = tb[as.character(x)]
  y[is.na(y)] = 0
  df1 = data.frame(type = 'Mapping-based', len = x, cnt = y)

  tb = table(nchar(tas$alt) - nchar(tas$ref))
  y = tb[as.character(x)]
  y[is.na(y)] = 0
  df2 = data.frame(type = 'Assembly-based', len = x, cnt = y)

  dos = cbind(rbind(df1, df2), org = org)
  do = rbind(do, dos)
}
#  do$cnt = log(do$cnt)

fo = sprintf("%s/compstat/comp_idm_size.pdf", Sys.getenv("misc3"))
p = ggplot(do) +
  geom_bar(mapping = aes(x = len, y = cnt, fill = type), 
    stat = 'identity', position = 'dodge', 
    geom_params=list(width = 0.8, alpha = 0.8)) +
  scale_fill_brewer(palette='Set1', name = '', guide = guide_legend(nrow = 1, byrow = T, label.position = "right", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  scale_x_continuous(name = 'Indel Size (bp)') +
  scale_y_continuous(name = '# events (log)') +
  facet_wrap( ~ org, nrow = 2) +
  theme_bw() +
  theme(legend.position = "top", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(0,1,1,0), "lines")) +
  theme(axis.text.x = element_text(size = 8, colour = "grey", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0))
ggsave(p, filename = fo, width = 6, height = 6)

##### validate SNP calling set
org = "HM340"
  dira = sprintf("%s/%s_HM101/23_blat/31.9", Sys.getenv('misc3'), org)
  fa = file.path(dira, "vnt.tbl")
  ta = read.table(fa, sep = "\t", header = F, as.is = T)
  colnames(ta) = c("chr", "pos", "ref", "alt", "score")

  dirm = file.path(Sys.getenv("misc3"), "hapmap", "12_ncgr", "44_tbl")
  fm = sprintf("%s/%s.tbl", dirm, org)
  tm = read.table(fm, sep = "\t", header = F, as.is = T)
  colnames(tm) = c("chr", "pos", "ref", "alt", "score")

  tms = tm[tm$chr %in% chrs & nchar(tm$ref) == 1 & nchar(tm$alt) == 1, ]
  tas = ta[ta$chr %in% chrs & nchar(ta$ref) == 1 & nchar(ta$alt) == 1, ]

dir = file.path(Sys.getenv("misc3"), "seqvalidation")
fg = file.path(dir, "16.snp")
tg = read.table(fg, sep = "\t", header = F, as.is = T)[,1:4]
colnames(tg) = c("chr", "pos", "ref", "alt")
tg = tg[order(tg$chr, tg$pos),]
tg = unique(tg)
tg = cbind(tg, id = paste(tg$chr, tg$pos, sep = ":"))
dup_ids = tg$id[duplicated(tg$id)]
tg = tg[!tg$id %in% dup_ids,1:4]


  tc = merge(tg, tas, by = c('chr', 'pos'), all.x = T)
  t_ovl = tc[!is.na(tc$alt.x) & !is.na(tc$alt.y),]
  t_mis = tc[!is.na(tc$alt.x) & is.na(tc$alt.y),]
  n_con = sum(t_ovl$alt.x == t_ovl$alt.y)
  n_dis = sum(t_ovl$alt.x != t_ovl$alt.y)
  n_con
  n_dis
  nrow(t_mis)

  tc = merge(tg, tms, by = c('chr', 'pos'), all.x = T)
  t_ovl = tc[!is.na(tc$alt.x) & !is.na(tc$alt.y),]
  t_mis = tc[!is.na(tc$alt.x) & is.na(tc$alt.y),]
  n_con = sum(t_ovl$alt.x == t_ovl$alt.y)
  n_dis = sum(t_ovl$alt.x != t_ovl$alt.y)
  n_con
  n_dis
  nrow(t_mis)
