require(plyr)
require(dplyr)
require(ggplot2)

dirw = '/home/springer/zhoux379/data/genome/Mo17/trinity'

fl = file.path(dirw, "02.len.tsv")
tl = read.table(fl, sep = "\t", as.is = T, header = F)
colnames(tl) = c('gid','size')
nrow(tl)


fb = file.path(dirw, "06.tsv")
tb = read.table(fb, sep = "\t", as.is = T, header = T)

tb2 = tb[tb$alnLen / tb$qSize >= 0.8 & tb$ident >= 0.8,]
length(unique(tb2$qId))

tlu = tl[!tl$gid %in% tb2$qId,]
nrow(tlu)
summary(tlu$size)