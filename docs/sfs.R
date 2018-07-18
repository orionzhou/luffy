#configurations on mthap
dir = '/export/lab/analysis/Selection/cds_HM102_20k_40_56fasta'
#setwd(dir)
f_id = file.path(dir, '02_ids.tbl')
dirI = file.path(dir, 'SILDIR')
#dirI = file.path(dir, 'REPDIR')
file_suffix = "all_silent"
#file_suffix = "replacement"
f_out = paste("sfs_", file_suffix, ".tbl", sep="")
dirE = '/export/lab/analysis/peng/Scripts/pub_c/bin'

#configurations on MSI
#dirI = '/project/youngn/zhoup/Data/in/SFStest'
#dirE = '/project/youngn/zhoup/Scripts/pub_c/bin'
#f_id = '/project/youngn/zhoup/Data/misc1/seq06/opt12/02_ids.tbl'

ids = read.table(f_id, sep="\t", header=T, as.is=T)$id
#ids = c("Medtr1g101010.1", "Medtr5g005010.1", "Medtr5g005030.1", "Medtr6g051860.1", "Medtr7g032360.1")
cutoff_missing = 55 # only consider positions where less than 'cutoff_missing' accessions have 'N's

cat(paste("id", "pos", "n_states", "anc", "der", "n_N", "n_anc", "n_der", "freq_der\n", sep="\t"), file=f_out, append=FALSE)
i = 0
for (id in ids) {
  fi = file.path(dirI, paste(id, ".fas.", file_suffix, sep=''))
  fo = "stat.tbl"
  if(!file.exists(fi)) next
  cmd = paste(file.path(dirE, "sspStat"), "-i", fi, "-o", fo, sep=' ')
  system(cmd)
  s01 = read.table(fo, header=T, sep="\t", as.is=T, colClasses=c(anc="character", der="character"))
  if(nrow(s01) == 0) next
  s02 = cbind(id=id, s01, freq_der = s01$n_der / (s01$n_anc+s01$n_der))
  s03 = s02[s02$n_states == 2 & s02$n_N < cutoff_missing,]
  write.table(s03, file=f_out, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  i = i+1
}
system("rm stat.tbl")


