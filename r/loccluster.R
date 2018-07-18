require(plyr)
require(ggplot2)
source("Location.R")

dirw = file.path(Sys.getenv("genome"), "HM056")
fg = file.path(dirw, '51.gtb')
tg = read.table(fg, header = T, as.is = T, sep = "\t", quote = "")

tc = tg[tg$cat2 == "NCR",]
tn = tg[tg$cat2 %in% c("CC-NBS-LRR", "TIR-NBS-LRR"),]

chrs = unique(tg$chr)
chrs = chrs[! chrs %in% c("chrU")]
dchr = data.frame(chrn = 1:length(chrs), chr = chrs, stringsAsFactors = F)

tx = tc
tx = tx[tx$chr %in% chrs,]

tx = merge(tx[,c('id','chr','beg','end')], dchr, by = 'chr')
tx = within(tx, {pos = chrn * 100000000 + (beg + end) / 2})
tx = tx[order(tx$pos), ]

dcl = locCluster(tx$pos, 100000)
to1 = ddply(dcl, .(cluster), summarise, cnt = length(id))
to2 = table(to1$cnt)
to3 = data.frame(clu_size = as.numeric(names(to2)), clu_num = as.numeric(to2), stringsAsFactors = F)
to3 = within(to3, {tot_num = clu_size * clu_num})
sum(to3$tot_num[to3$clu_size > 2]) / sum(to3$tot_num)

