dir = file.path(DIR_Misc2, "crp")
d03 = file.path(dir, "03_hits")
d04 = file.path(dir, "04_hits_picked")
d05 = file.path(dir, "05_add_col")
d08 = file.path(dir, "08_curation")
d09 = file.path(dir, "09_assembled")
d11 = file.path(dir, "11_model")

f04_20 = file.path(d04, "20_final.tbl")

f05_02 = file.path(d05, "02_status.tbl")
f05_13 = file.path(d05, "13_db_model_3.tbl")

f08_01 = file.path(d08, "01.tbl")
f08_02 = file.path(d08, "02.tbl")
f08_03 = file.path(d08, "03_ks.tbl")
f08_04 = file.path(d08, "04_pz.tbl")

f09_01 = file.path(d09, "01.tbl")
f09_02 = file.path(d09, "02.tbl")
f09_03 = file.path(d09, "03.tbl")
f09_04 = file.path(d09, "04_cleaned.tbl")

assembly1 = function(f09_01) {
  th = read.table(f04_20, sep="\t", as.is=T, quote="", header=T)
  tm = read.table(f05_13, sep="\t", as.is=T, quote="", header=T)
  tc = read.table(f08_01, sep="\t", as.is=T, quote="", header=T)
  ts = read.table(f05_02, sep="\t", as.is=T, quote="", header=T)

  tm.1 = tm[,c('chr','location','id','tag')]
  colnames(tm.1) = c('chr','location','hit_id','hit_tag')

  t.2 = merge(th, tm.1, by=c('chr','location'))
  t.3 = t.2

  tc.1 = tc[, c("id1", "comment_new")]
  idxs = is.na(tc.1$comment_new)
  tc.1$comment_new[idxs] = tc$comment[idxs]
  colnames(tc.1)[2] = "comment"
  t.4 = merge(t.3, tc.1, by='id1')

  t.5 = cbind(t.4, comment_new="", comment_more="")

  ts.1 = ts[,c('id1','status')]
  colnames(ts.1)[2] = "comment_old"
  t.6 = merge(t.5, ts.1, by='id1')
  
  t.7 = t.6[order(t.6$id1),]
  write.table(t.7, f09_01, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, na="")
}
assembly2 = function(f09_03) {
  th = read.table(f04_20, sep="\t", as.is=T, quote="", header=T)
  tm = read.table(f05_13, sep="\t", as.is=T, quote="", header=T)

  tm.1 = tm[,c('chr','location','id','tag')]
  colnames(tm.1) = c('chr','location','hit_id','hit_tag')

  t.2 = merge(th, tm.1, by=c('chr','location'))

  f_mtb = file.path(DIR_Misc2, "mapping/16_mt_35_crp/12_rmdup.mtb")
  m01 = read.table(f_mtb, sep="\t", quote="", header=T, as.is=T)
  m01.1 = m01[,-1]

  t.5 = merge(t.2, m01.1, by.x='id1', by.y='qId', all.x=T)

  t.9 = t.5[order(t.5$id1),]
  write.table(t.9, f09_03, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, na="")
}
assembly3 = function(f09_04) {
  th = read.table(f04_20, sep="\t", header=T, as.is=T)
  tc = read.table(f08_04, sep="\t", header=T, as.is=T, quote="")[,c(1,6:8,11:12)] 
  t.2 = merge(th, tc, by='id1')

  t.5 = t.2[t.2$family < "CRP1600",]
  t.6 = t.5[t.5$comment3 != 1,]
  cat(sprintf("%d false positive\n", nrow(t.5)-nrow(t.6)))
  table(t.6$comment3)
  t.8=t.6[order(t.6$id1),c(1:5,11:13)]
  write.table(t.8, f09_04, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, na="")
}


dir = file.path(DIR_Misc2, "crp/03_hits")
fo = file.path(dir, "12_picked_201112.tbl")
fn = file.path(dir, "12_picked.tbl")
to = read.table(fo, header=T, sep="\t", quote="", as.is=T)
tn = read.table(fn, header=T, sep="\t", quote="", as.is=T)

tr = cbind(tn[,c('id1', 'chr', 'beg', 'end', 'strand')], n_id2=NA, id2=NA)
for (i in 1:nrow(tn)) {
  chr = tn$chr[i]
  beg = tn$beg[i]
  end = tn$end[i]
  str = tn$strand[i]
  len_off = 30
  to.1 = to[to$chr==chr & to$strand==str & abs(to$beg-beg) < len_off & abs(to$end-end) < len_off, ]
  tr$n_id2[i] = nrow(to.1)
  tr$id2[i] = to.1$id1[1]
}
table(tr$n_id2)

