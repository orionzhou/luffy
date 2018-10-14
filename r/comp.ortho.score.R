require(rtracklayer)
require(plyr)
require(seqinr)
source("Align.R")
source("comp.fun.R")

#args <- commandArgs(trailingOnly = TRUE)
#print(args)
#qname <- ifelse(is.na(args[1]), "HM004", as.character(args[1]))

##### cluster initialization
cl = makeCluster(detectCores())

cluster_fun <- function() {
    require(Biostrings)
    require(plyr)
}
clusterCall(cl, cluster_fun)

##### compute ortholog identity using pairwise alignment
dirw = file.path(Sys.getenv("misc3"), "comp.ortho.hm")

f_tfas = file.path(Sys.getenv("genome"), tname, "51.fas")
tfas <- read.fasta(f_tfas, seqtype = "AA", as.string = T, set.attributes = F)
tids = names(tfas)

qnames = qnames_15
for (qname in qnames) {
  #qname = "HM004"
  f_qfas = file.path(Sys.getenv("genome"), qname, "51.fas")
  qfas <- read.fasta(f_qfas, seqtype = "AA", as.string = T, set.attributes = F)
  qids = names(qfas)

  fi = sprintf("%s/01_syn_ortho/%s.tbl", dirw, qname)
  ti = read.table(fi, sep = "\t", header = T, as.is = T)[,1:5]
  ti = ti[ti$qid != '' & ti$qid %in% qids & ti$tid %in% tids,]

  tm = cbind(ti, tseq = as.character(tfas[ti$tid]), 
    qseq = as.character(qfas[ti$qid]), stringsAsFactors = F)

  ptm <- proc.time()
  y = parApply(cl, tm, 1, aa_pw_dist)

  cat(qname, proc.time() - ptm, "\n")
  
  to = cbind(ti, t(y))
  to$qlen = as.integer(to$qlen / 3)
  to$tlen = as.integer(to$tlen / 3)
  fo = sprintf("%s/11_score/%s.tbl", dirw, qname)
  write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
}

##### stop cluster
stopCluster(cl)