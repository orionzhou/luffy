require(plyr)
require(dplyr)
require(tidyr)

dirw = file.path(Sys.getenv('misc2'), 'wheatped')


fi = file.path(dirw, '12.inb.coef.tsv')
ti = read.table(fi, header = F, sep = "\t", as.is = T)

sids = unique(ti$V1[!startsWith(ti$V1, "ped")])
nsid = length(sids)
cop = matrix(rep(NA, nsid*nsid), nrow=nsid)
colnames(cop) = sids
rownames(cop) = sids

# cut -f1-3 24.pair.cop/* > 24.pair.cop.tsv
f1 = sprintf("%s/24.pair.cop.tsv", dirw)
t1 = read.table(f1, header = F, sep = "\t", as.is = T)
t1 = t1[order(t1$V2, t1$V2),]
t1 = t1[!duplicated(t1[,1:2]),]
t2 = spread(t1, V2, V3)
t3 = t2[,-1]
rownames(t3) = t2$V1

missingcol = sids[!sids %in% colnames(t3)]
missingcol
t3 = cbind(t3, NA)
colnames(t3)[ncol(t3)] = missingcol

missingrow = sids[!sids %in% rownames(t3)]
missingrow
t3 = rbind(t3, NA)
rownames(t3)[nrow(t3)] = missingrow

t4 = t3[sids, sids]
tmp = t4[upper.tri(t4)]
t4 = t(t4)
t4[upper.tri(t4)] = tmp

diag(t4) = 1
sum(is.na(t4))

fo = file.path(dirw, "26.cop.tsv")
write.table(t4, fo, sep = "\t", row.names = T, col.names = T, quote = F)
