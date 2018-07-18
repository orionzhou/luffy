require(Gviz)
require(rtracklayer)
source("comp.fun.R")
source("gviz.fun.R")

dirw = file.path(Sys.getenv("misc1"), 'linhao')

# blat $genome/HM340.AC/11_genome.fas 01.fas 11.ac.psl
# blat $genome/HM340.FN/11_genome.fas 01.fas 12.fn.psl
# psl2gal.pl -i 11.ac.psl -o 11.ac.gal
# psl2gal.pl -i 12.fn.psl -o 12.fn.gal

fi = file.path(dirw, "11.ac.gal")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti[order(ti$ali, decreasing = T),][1:20,c(2:11,13:15)]

fi = file.path(dirw, "12.fn.gal")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti[order(ti$ali, decreasing = T),][1:20,c(2:11,13:15)]

chr7:13900000-14600000 and one 100kb gap at chr7:14435000-14535000)

scf016:1820965-1829330
scf016:1820001-1830000