library(ggplot2)
library(plyr)
library(data.table)
library(RColorBrewer)

ff = file.path(DIR_Misc2, "crp.annotation/family_info.tbl")
tf = read.table(ff, header=TRUE, sep="\t", as.is=T)
fams = unique(tf$cat)

famstrs = c()
for (fam in fams) {
  idxs= which(tf$cat == fam)
  idx_idxs_b = which( idxs != c(-1, idxs[-length(idxs)]+1) )
  idx_idxs_e = c(idx_idxs_b[-1]-1, length(idxs))

  fameles = c()
  for (j in 1:length(idx_idxs_b)) {
    famele = ifelse( idx_idxs_b[j] == idx_idxs_e[j], tf$id[idxs[idx_idxs_b[j]]], 
      sprintf("%s-%s", tf$id[idxs[idx_idxs_b[j]]], tf$id[idxs[idx_idxs_e[j]]]) )
    fameles = c(fameles, famele)
  }
  famstr = paste(fameles, collapse=",")
  famstrs = c(famstrs, famstr)
}
fami = data.frame(cat=fams, ids=famstrs)

diro = file.path(DIR_Drop, "Docs/research/MS_SPADA/supplement")

org = "Mtruncatula_4.0"
orgs = c("Athaliana", "Mtruncatula_3.5")
pres = c("crp", "defl")
pres = c("crp", "crp")

dc=c()
for (i in 1:length(orgs)) {
  org = orgs[i]
  fg = sprintf("%s/spada.%s.%s/31_model_evaluation/61_final.tbl", DIR_Misc4, pres[i],  org)
  tg = read.table(fg, header=TRUE, sep="\t", quote="", as.is=T)[,c("id","family")]
  tg2 = merge(tg, tf, by.x='family', by.y='id')
  tg3 = ddply(tg2, .(cat), "nrow")
  colnames(tg3)[2] = org
  if(i==1) {
    dc = tg3
  } else {
    dc = merge(dc, tg3, by.x='cat', all=T)
  }
}
dc = merge(fami, dc, by="cat", all=T)
dc[is.na(dc)] <- 0
colnames(dc)[1] = c('family')
dc = dc[order(dc$ids),]
apply(dc[,orgs], 2, sum)
write.table(dc, file.path(diro, 'crp.family.number.tbl'), col.names=T, row.names=F, quote=F, sep='\t')



org = orgs[2]
fg = sprintf("%s/spada.crp.%s.simple/31_model_evaluation/61_final.tbl", DIR_Misc4, org)
t1 = read.table(fg, header=TRUE, sep="\t", quote="", as.is=T)
colnames(t1)[which(colnames(t1)=="start")] = "beg"
t2 = t1[tolower(substr(t1$chr,1,3))=="chr" & substr(t1$chr,4,5)<=9 & substr(t1$chr,4,5)>=0, -which(colnames(t1)=="sequence")]

tli = t2
locs = as.numeric(substr(tli$chr, 4, 5)) * 1000000000 + (tli$beg + tli$end) / 2
names(locs) = tli$id
cl = locCluster(locs, 50000)
tlo = merge(tli, cl, by="id")

td = merge(tlo, tf, by.x="family", by.y="id")

pal <- colorRampPalette(brewer.pal(9,"Set1"))
p <- ggplot(data=td) +
  geom_point(mapping=aes(x=(beg+end)/2/1000000, y=cluster_y, colour=factor(cat, levels=fams)), shape=18, size=0.8) +
  scale_colour_manual(values=pal(length(unique(td$cat)))) +
  scale_x_continuous(name='Chr Position (Mbp)', expand=c(0.01, 0)) +
  scale_y_continuous(name='', expand=c(0.04, 0)) +
  facet_grid(chr ~ .) + 
  theme(legend.position='right', legend.title=element_blank()) +
  theme(axis.text.x=element_text(size=8, angle=0)) +
  theme(axis.text.y=element_blank(), axis.ticks=element_blank())
ggsave(p, filename=sprintf("%s/%s.pdf", diro, org), width=8, heigh=5)

#write.table(g21[order(g21$pos),], file.path(dir, "02_cluster.tbl"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

org = "Athaliana"
org = "Osativa"
fs1 = file.path(DIR_Misc2, "crp.ssp", org, "11_stat_una.tbl")
ts1 = read.table(fs1, header=T, sep="\t", as.is=T)
table(ts1$tag)
fs2 = file.path(DIR_Misc2, "crp.ssp", org, "11_stat_spa.tbl")
ts2 = read.table(fs2, header=T, sep="\t", as.is=T)
table(ts2$tag)
fs3 = file.path(DIR_Misc2, "crp.ssp", org, "11_stat_cur.tbl")
ts3 = read.table(fs3, header=T, sep="\t", as.is=T)
table(ts3$tag)

org = "Athaliana"
org = "Mtruncatula_3.5"
org = "Osativa"
fi = sprintf("%s/spada.crp.%s.simple/33_stat_spa.tbl", DIR_Misc4, org)
ti = read.table(fi, sep="\t", header=T, as.is=T)
table(ti$tag)

orgs = c("Mtruncatula_4.0", "HM056", "HM340")
for (org in orgs) {
  fi = sprintf("%s/spada.crp.%s/31_model_evaluation/61_final.tbl", DIR_Misc4, org)
  ti = read.table(fi, sep="\t", header=T, as.is=T)[,1:6]
  n_crp = nrow(ti)
  n_ncr = sum(ti$family >= "CRP1130" & ti$family <= "CRP1530")
  cat(org, "CRP:", n_crp, "NCR:", n_ncr, "\n", sep=" ")
}


org = "HM101"
fi = sprintf("/home/youngn/zhoup/Data/misc4/spada.crp.%s/61_final.tbl", org)
ti = read.table(fi, sep="\t", header=T, as.is=T)[,c(1:5)]
ts = ti[ti$family <= 'CRP1600',]

