#configurations on mthap
dirI = '/export/lab/analysis/Selection/SFStest' 
dirE = '/export/lab/analysis/peng/Scripts/pub_c/bin'
f_id = '/export/lab/analysis/peng/Data/in/02_ids.tbl'

#configurations on MSI
dirI = '/project/youngn/zhoup/Data/in/SFStest'
dirE = '/project/youngn/zhoup/Scripts/pub_c/bin'
f_id = '/project/youngn/zhoup/Data/misc1/seq06/22_seq/02_ids.tbl'

ids = read.table(f_id, sep="\t", header=T, as.is=T)$id
intervals = seq(0,1,length.out=20)  # initialize 20 bins
cutoff_missing = 16 # only consider positions where less than 'cutoff_missing' accessions have 'N's

d01 = matrix(NA, nrow=0, ncol=length(intervals)-1)
colnames(d01) = names(table(cut(0,breaks=intervals)))
i = 0
for (id in ids) {
  fi = file.path(dirI, paste(id, ".fas.all_silent", sep=''))
  fo = "stat.tbl"
  
  if(!file.exists(fi)) next
  i = i+1
  cmd = paste(file.path(dirE, "sspStat"), "-i", fi, "-o", fo, sep=' ')
  system(cmd)
  s01 = read.table(fo, header=T, sep="\t", as.is=T)
  s02 = cbind(s01, freq_der = s01$n_der / (s01$n_anc+s01$n_der))
  s03 = s02[s02$n_states == 2 & s02$n_N < cutoff_missing,]

  t01 = table(cut(s03$freq_der, breaks=intervals))
  d01 = rbind(d01, t01)
  rownames(d01)[i] = id
}
write.table(d01, 'sfs_silent.tbl', sep="\t", quote=F, row.names=T, col.names=T)

d02 = matrix(NA, nrow=0, ncol=length(intervals)-1)
colnames(d02) = names(table(cut(0,breaks=intervals)))
i = 0
for (id in ids) {
  fi = file.path(dirI, paste(id, ".fas.replacement", sep=''))
  fo = "stat.tbl"

  if(!file.exists(fi)) next
  i = i+1
  cmd = paste(file.path(dirE, "sspStat"), "-i", fi, "-o", fo, sep=' ')
  system(cmd)
  s01 = read.table(fo, header=T, sep="\t", as.is=T)
  s02 = cbind(s01, freq_der = s01$n_der / (s01$n_anc+s01$n_der))
  s03 = s02[s02$n_states == 2 & s02$n_N < cutoff_missing,]

  t01 = table(cut(s03$freq_der, breaks=intervals))
  d02 = rbind(d02, t01)
  rownames(d02)[i] = id
}
write.table(d02, 'sfs_replace.tbl', sep="\t", quote=F, row.names=T, col.names=T)

system("rm stat.tbl")


