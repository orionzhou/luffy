source("Location.R")

  orgs = c(
    "HM058", "HM056", "HM125", "HM129", "HM034", 
    "HM095", "HM060", "HM185", "HM004", "HM050", 
    "HM023", "HM010", "HM022", "HM340", "HM324"
  )

for (org in orgs) {
#org = "HM101"

dirw = file.path(Sys.getenv("genome"), org)
dir.create(dirw, showWarnings = F)
setwd(dirw)

f_fas = "11_genome.fas"
stopifnot(file.exists(f_fas))

fis = c("augustus/41.gtb", "42.nbs/11.gtb", "44.rlk/11.gtb")
fos = c("41.gtb", "42.nbs.gtb", "44.rlk.gtb")

for (i in 1:length(fis)) {
  fi = fis[i]; fo = fos[i]
  stopifnot(file.exists(fi))
  if(file.exists(fo)) file.remove(fo)
  file.symlink(fi, fo)
}

f_crp = sprintf("%s/spada.crp.%s/61_final.gtb", Sys.getenv('misc4'), org)
stopifnot(file.exists(f_crp))
tc = read.table(f_crp, header = T, sep = "\t", as.is = T, quote = "")[,1:18]
tc$cat1 = "mRNA"
tc$note = ""
tc$cat2[tc$cat3 >= "CRP0000" & tc$cat3 <= "CRP1030"] = "CRP0000-1030"
tc$cat2[tc$cat3 >= "CRP1040" & tc$cat3 <= "CRP1530"] = "NCR"
tc$cat2[tc$cat3 >= "CRP1600" & tc$cat3 <= "CRP6250"] = "CRP1600-6250"
write.table(tc, "43.crp.gtb", sep = "\t", row.names = F, col.names = T, quote = F, na = '')

system("gtb.merge.pl -a 41.gtb -b 42.nbs.gtb -o 49.1.gtb")
system("gtb.merge.pl -a 49.1.gtb -b 43.crp.gtb -o 49.2.gtb")
system("gtb.merge.pl -a 49.2.gtb -b 44.rlk.gtb -o 49.gtb")
system("gtb.dedup.pl -i 49.gtb -o 50.1.dedup.gtb")

### refine TE annotation
fr = "12.rm.tbl"
stopifnot(file.exists(fr))
tr = read.table(fr, sep = "\t", header = F, as.is = T)
idxs = grep("^(DNA)|(LINE)|(LTR)|(SINE)|(RC)", tr$V9)
grr = with(tr[idxs,], GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
grr = reduce(grr)

tg = read.table("50.1.dedup.gtb", sep = "\t", header = T, as.is = T, quote = "")
idxs = which(tg$cat2 == 'Unknown')
gr = with(tg[idxs,], GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

bps = intersect_basepair(gr, grr)
pct = bps / width(gr)
idxs_te = idxs[pct >= 0.7]
tg$cat2[idxs_te] = 'TE'
tg$note[idxs_te] = 'RepeatMasker'
write.table(tg, "50.2.te.gtb", sep = "\t", row.names = F, col.names = T, quote = F, na = '')

### routine
system("gtb.pickalt.pl -i 50.2.te.gtb -o 50.3.pickalt.gtb")

if(org != "HM101") {
  f_ctm = sprintf("%s/%s_HM101/41_novseq/15.foreign.scf.txt", Sys.getenv("misc3"), org)
  stopifnot(file.exists(f_ctm))
  system(sprintf("gtb.fill.len.pl -i 50.3.pickalt.gtb | gtb.filter.pl -l 30 -c %s -o 51.gtb", f_ctm))
} else {
  system("gtb.fill.len.pl -i 50.3.pickalt.gtb | gtb.filter.pl -l 30 -o 51.gtb")
}

tg = read.table("51.gtb", sep = "\t", header = T, as.is = T, quote = "")[,c(1:6,16:18)]
table(tg$cat2)

}


for (org in orgs) {
dirw = file.path(Sys.getenv("genome"), org)
dir.create(dirw, showWarnings = F)
setwd(dirw)

system("gtb2gff.pl -i 51.gtb -o 51.gff")
system("gtb.idx.pl -i 51.gtb -s 15.sizes")
system("gtb2tbl.pl -i 51.gtb -o 51.tbl")
system("gtb2fas.pl -i 51.gtb -d 11_genome.fas -o 51.fas")
}