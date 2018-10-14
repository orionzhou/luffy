require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)

dirw = file.path(Sys.getenv('misc1'), 'expdesign')

fi = file.path(dirw, 'tissues.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

genos = c("B73", "Mo17", "B73xMo17", "Mo17xB73")

to = ddply(ti, .(TissueID), expand <- function(x, genos) {
	data.frame(genotype=genos[1:x[1,'Genotypes']])
}, genos = genos)
to$TissueID = factor(to$TissueID, levels = ti$TissueID)
to$genotype = factor(to$genotype, levels = genos)
to = to[order(to$TissueID, to$genotype), ]

to = data.frame(Tissue = rep(to$TissueID, each = 3), 
	Genotype = rep(to$genotype, each = 3),
	Replicate = rep(c(1:3), nrow(to)))

fo = file.path(dirw, "05.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
