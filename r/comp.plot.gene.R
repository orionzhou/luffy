require(rtracklayer)
require(xlsx)
source("comp.fun.R")
source("comp.plot.fun.R")

dirw = file.path(Sys.getenv("misc3"), 'comp.plot.gene')
fl = file.path(dirw, 'loci.xlsx')
tl = read.xlsx(fl, sheetIndex = 1, header = T, stringsAsFactors = F)

fid = file.path(Sys.getenv("misc3"), "comp.ortho.hm", "01.ids.tbl")
tid = read.table(fid, header = T, sep = "\t", as.is = T)

fds = file.path(Sys.getenv("misc3"), "comp.ortho.hm", "12.score.tbl")
tds = read.table(fds, header = T, sep = "\t", as.is = T)

fg = file.path(Sys.getenv("genome"), "HM101", "51.gtb")
tg = read.table(fg, header = T, sep = "\t", as.is = T)

fr = file.path(Sys.getenv("misc3"), "comp.og/05.clu/32.tbl")
tr = read.table(fr, header = T, sep = "\t", as.is = T)
orgs = sapply(strsplit(tr$id, split="-"), "[", 1)
gids = sapply(strsplit(tr$id, split="-"), "[", 2)
tr = cbind(tr, org = orgs, gid = gids)

##### 
source("comp.plot.fun.R")

i = 1
chr = tl$chr[i]; beg = tl$beg[i]; end = tl$end[i]
gidxs = which(tg$chr == chr & tg$beg >= beg & tg$end <= end)



fn = sprintf("%s/fig%03d.pdf", dirw, i)
CairoPDF(file = fn, width = 7, height = res$ht/72, bg = 'transparent')
grid.newpage()
grid.draw(res$grobs)
dev.off()


