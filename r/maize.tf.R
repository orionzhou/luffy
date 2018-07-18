require(tidyverse)
require(GenomicRanges)

dirw = file.path(Sys.getenv("genome"), "Zmays_v4/TF")

### 
fg = file.path(dirw, "../51.gtb")
tg = read_tsv(fg, col_types = "ccciiccccccccccccc") %>%
    transmute(tid = id, gid = par, chrom = chr)

fm = file.path(dirw, "../gene_mapping/maize.v3TOv4.geneIDhistory.txt")
tm = read_tsv(fm, col_names = F) %>%
    transmute(ogid = X1, gid = X2, change = X3, method = X4, type = X5)

Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

fi = file.path(dirw, 'Zma_TF_list')
ti = read_tsv(fi, col_names = T) %>%
    transmute(ogid = Gene_ID, otid = TF_ID, fam = Family) %>%
    filter(ogid != '') %>%
    group_by(ogid) %>%
    summarise(fam = Mode(fam))

tf = ti %>%
    inner_join(tm, by = 'ogid') %>%
    filter(gid %in% tg$gid) %>%
    #filter(type == '1-to-1') %>%
    select(gid, fam)
tfs = tf %>% distinct(fam) %>% arrange(fam) %>%
    mutate(fid = sprintf("tf%04d", 1:length(fam)))
tf = tf %>% inner_join(tfs, by = 'fam') %>%
    select(fid, gid, fam) %>%
    arrange(fid, gid)

tf %>% count(fid)
fo = file.path(dirw, "10.tsv")
write_tsv(tf, fo)


### find house keeping genes
dirw = file.path(Sys.getenv("genome"), "Zmays_v4/housekeeping")

fi = file.path(dirw, '01.tsv')
ti = read.table(fi, header = F, sep = "\t", as.is = T)
colnames(ti) = c("ogid", 'name')

ti2 = merge(ti, tm, by = 'ogid')
sum(ti2$ngid %in% tg$par)

fo = file.path(dirw, "11.tsv")
write.table(ti2, fo, sep = "\t", row.names = F, col.names = T, quote = F)
