require(tidyverse)
require(GenomicRanges)

fg = file.path(dirw, "../51.gtb")
tg = read_tsv(fg, col_types = 'ccciiccccccccccccc')
gids = unique(tg$par)

# use interpro
dirg = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro'

fg = file.path(dirg, "../51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,1:2]
colnames(tg) = c("tid", "gid")

fm = file.path(dirg, "../57.longest.tsv")
tm = read.table(fm, sep = "\t", header = F, as.is = T)
colnames(tm) = c("gid", 'tid')

ft = file.path(dirg, "07.tsv")
tt = read.table(ft, sep = "\t", header = F, as.is = T)
colnames(tt) = c("tid", "gos")

tt2 = merge(tg, tt, by = 'tid', all.x = T)
stopifnot(nrow(tg) == nrow(tt2))
tt2 = tt2[order(tt2$gid, tt2$tid),]

fo1 = file.path(dirg, "10_mrna.tsv")
write.table(tt2[,c(1,3)], fo1, sep = "\t", row.names = F, col.names = F, quote = F, na = '')

tm2 = merge(tm, tt, by = 'tid', all.x = T)
stopifnot(nrow(tm2) == nrow(tm))
tm2 = tm2[order(tm2$gid, tm2$tid),]
fo2 = file.path(dirg, "11_gene.tsv")
write.table(tm2[,c(2,3)], fo2, sep = "\t", row.names = F, col.names = F, quote = F, na = '')

# process V4 GO
dirw = '/home/springer/zhoux379/data/genome/Zmays_v4/GO'
fgo = '/home/springer/zhoux379/data/db/go/go-basic.tsv'
tgo = read_tsv(fgo) %>%
    transmute(goid = goid, level = level, goname = name)

to = tibble()
ctags = c("Interproscan5", "arabidopsis", "argot2.5", "fanngo", 
          "pannzer", "uniprot.plants", "aggregate")
for (ctag in ctags) {
    fi = sprintf("%s/maize.B73.AGPv4.%s.gaf", dirw, ctag)
    ti = read_tsv(fi, skip = 1)[,c(2,5,7,9)]
    colnames(ti) = c("gid", "goid", "evidence", "gotype")
    ti = ti %>% mutate(ctag = ctag) %>%
        select(ctag, gid, goid, evidence, gotype)
    to = to %>% bind_rows(ti)
}

goids = unique(to$goid)
length(goids)
sum(goids %in% tgo$goid)
to2 = to %>% inner_join(tgo, by = 'goid')
to2 %>% count(ctag)
to2 %>% distinct(ctag, gid) %>% count(ctag)

fo = file.path(dirw, "09.tsv")
write_tsv(to2, fo)

# add corncyc / TF
dirw = '/home/springer/zhoux379/data/genome/Zmays_v4/GO'
fi = file.path(dirw, "09.tsv")
ti = read_tsv(fi)

fp = file.path(dirw, '../corncyc/10.tsv')
tp = read_tsv(fp) %>%
    transmute(ctag = 'corncyc', gid = gid, goid = pid,
              evidence = '', gotype = '', level = '', goname = pathway)

ff = file.path(dirw, '../TF/10.tsv')
tf = read_tsv(ff) %>%
    transmute(ctag = 'tfdb', gid = gid, goid = fid,
              evidence = '', gotype = '', level = '', goname = fam)

to = rbind(ti, tp, tf)

to %>% count(ctag)
to %>% distinct(ctag, gid) %>% count(ctag)
fo = file.path(dirw, "10.tsv")
write_tsv(to, fo)

