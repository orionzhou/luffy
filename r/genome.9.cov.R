require(GenomicRanges)
require(dplyr)

dirg = '/home/springer/zhoux379/scratch/mo17vnt/51.realign'

fc = file.path(dirg, "Mo17.cov.tsv")
tc = read.table(fc, sep = "\t", header = F, as.is = T)
tc = tc[tc$V1 == 'genome',]

dat = tc[,2:3]
colnames(dat) = c('value', 'freq')
mu = sum(as.numeric(dat$freq * dat$value)) / sum(dat$freq)
sigma = with(dat, sqrt(sum(freq*(value-mu)^2)/(sum(freq)-1)))

cutoff = mu + 5 * sigma
sum(dat$value > mu + 5*sigma)

fp = file.path(dirg, "cov.pdf")
pdf(fp, height = 4, width = 5)
plot(dat$value, dat$freq/1000000, type = 'l', xlim = c(0, cutoff), xlab = 'Depth of Coverage', ylab = 'Million Bases')
dev.off()