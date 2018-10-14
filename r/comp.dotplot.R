########## load genome1 datasets
dat1.name = "hm056"
dat1.name = "hm340"

dat1.dir = file.path("/home/youngn/zhoup/Data/misc3", dat1.name)

# sequence lengths
f_len = file.path(dat1.dir, "11_seqlen.tbl")
dat1.seqlen = read.table(f_len, header=TRUE, sep="\t", as.is=T)

# gap locations
f_gap = file.path(dat1.dir, "12_gaploc.tbl")
dat1.gap = read.table(f_gap, header=TRUE, sep="\t", as.is=T)

######### load genome2 datasets
dat2.name = "hm101"

dat2.dir = "/home/youngn/zhoup/Data/genome/Mtruncatula_4.0"

# sequence lengths
f_len = file.path(dat2.dir, "15_seqlen.tbl")
dat2.seqlen = read.table(f_len, header=TRUE, sep="\t", as.is=T)

# gap locations
f_gap = file.path(dat2.dir, "16_gap_loc.tbl")
dat2.gap = read.table(f_gap, header=TRUE, sep="\t", as.is=T)


########## load (blast) comparison dataset
tc1 = read.table(file.path(dat1.dir, "23_blat/17_chain.gal"), header=TRUE, sep="\t", as.is=T)
tb1 = read.table(file.path(dat1.dir, "23_blat/18_block.gal"), header=TRUE, sep="\t", as.is=T)

tc = tc1
tb = tb1

qid.levels = unique(tc$qId[order(tc$hId, tc$hBeg)])
dat1.coord = data.frame(id=qid.levels, order=1:length(qid.levels))
dat1.coord = merge(dat1.coord, dat1.seqlen, by='id')
dat1.coord = dat1.coord[order(dat1.coord$order),-2]
dat1.coord = cbind(dat1.coord, beg=1, end=dat1.coord$length)
for (i in 2:nrow(dat1.coord)) {
  dat1.coord$beg[i] = dat1.coord$end[i-1] + 1
  dat1.coord$end[i] = dat1.coord$beg[i] + dat1.coord$length[i] - 1
}


dat2.coord = cbind(dat2.seqlen, beg=1, end=dat2.seqlen$length)
for (i in 2:nrow(dat2.coord)) {
  dat2.coord$beg[i] = dat2.coord$end[i-1] + 1
  dat2.coord$end[i] = dat2.coord$beg[i] + dat2.coord$length[i] - 1
}

alnplot = tb
colnames(dat1.coord) = c('qId', 'qLen.r', 'qBeg.r', 'qEnd.r')
alnplot = merge(alnplot, dat1.coord, by='qId')
colnames(dat2.coord) = c('hId', 'hLen.r', 'hBeg.r', 'hEnd.r')
alnplot = merge(alnplot, dat2.coord, by='hId')

app = alnplot[alnplot$qSrd == '+',]
apn = alnplot[alnplot$qSrd == '-',]
app = cbind(app, qBeg.a=app$qBeg-1+app$qBeg.r, qEnd.a=app$qEnd-1+app$qBeg.r, hBeg.a=app$hBeg-1+app$hBeg.r, hEnd.a=app$hEnd-1+app$hBeg.r)
apn = cbind(apn, qBeg.a=apn$qLen.r-apn$qEnd+apn$qBeg.r, qEnd.a=apn$qLen.r-apn$qBeg+apn$qBeg.r, hBeg.a=apn$hBeg-1+apn$hBeg.r, hEnd.a=apn$hEnd-1+apn$hBeg.r)

p <- ggplot(rbind(app,apn)) +
  geom_segment(mapping=aes(x=hBeg.a, xend=hEnd.a, y=qBeg.a, yend=qEnd.a), size=0.5) +
  layer(data=dat2.coord, geom='rect', mapping=aes(xmin=hBeg.r, xmax=hEnd.r, ymin=-3000000, ymax=0, fill=hId), geom_params=list(size=0)) +
  layer(data=dat2.coord, geom='text', mapping=aes(x=(hBeg.r+hEnd.r)/2, y=-4000000, label=hId), geom_params=list(hjust=0.5, vjust=1, angle=0, size=3)) +
  scale_fill_brewer(palette="Set1") +
  scale_x_continuous(name='HM101 (Mt4.0)') +
  scale_y_continuous(name=toupper(dat1.name)) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank())
ggsave(file.path(dat1.dir, "all.png"), p, width=8, height=8)