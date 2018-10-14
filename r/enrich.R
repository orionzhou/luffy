require(tidyverse)

#{{{ old interprosan GO output
if (FALSE) {
    fg = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/15.tsv'
    tg = read.table(fg, sep = "\t", as.is = T, header = F, quote = '')
    colnames(tg) = c("gid", "goid")

    fd = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/16.go.tsv'
    #td = read.table(fd, sep = "\t", header = T, as.is = T, quote = '')

    fi = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/gids.txt'
    ti = read.table(fi, as.is = T)
    gids.all = ti$V1
}
#}}}

fgo = '/home/springer/zhoux379/data/genome/B73/GO/10.tsv'
tgo = read_tsv(fgo)
unique(tgo$ctag)
gids_all = unique(tgo$gid)
tgo_ipr = tgo %>% filter(ctag == 'Interproscan5') %>% select(-ctag)
tgo_uni = tgo %>% filter(ctag == 'uniprot.plants') %>% select(-ctag)
tgo_ath = tgo %>% filter(ctag == 'arabidopsis') %>% select(-ctag)
tgo_arg = tgo %>% filter(ctag == 'argot2.5') %>% select(-ctag)

go_enrich <- function(gids, tg) {
    #{{{
    tgn = tg %>% distinct(goid, goname, gotype, level)
    tgs = tg %>% count(goid) %>% transmute(goid=goid, hitInPop = n)

    gids_all = tg %>% distinct(gid) %>% pull(gid)
    gids = unique(gids)
    gids = gids[gids %in% gids_all]
    sampleSize = length(gids)

    tz = tg %>% filter(gid %in% gids) %>%
        count(goid) %>% 
        transmute(goid = goid, hitInSample = n)

    tw = tz %>% inner_join(tgs, by = 'goid') %>% 
        mutate(sampleSize = length(gids), popSize = length(gids_all), 
            pval.raw = phyper(hitInSample-1, hitInPop, popSize-hitInPop, sampleSize, lower.tail = F),
            pval.adj = p.adjust(pval.raw, method = "BH")) %>%
        filter(pval.raw < 0.05) %>%
        arrange(pval.adj) %>%
        transmute(goid = goid, ratioInSample = sprintf("%d/%d", hitInSample, sampleSize),
            ratioInPop = sprintf("%d/%d", hitInPop, popSize),
            pval.raw = pval.raw, pval.adj = pval.adj)

    tw %>% left_join(tgn, by = 'goid') %>% arrange(pval.adj, pval.raw)
    #}}}
}

go_enrich_gosets <- function(gids, tgo. = tgo, pval.cutoff = 0.05,
                             srcs = c("uniprot.plants", "arabidopsis", "corncyc", "tfdb", "Interproscan5")) {
    #{{{
    to = tibble()
    for (src in srcs) {
        tgoc = tgo %>% filter(ctag == src) %>% select(-ctag)
        to1 = go_enrich(gids, tgoc) %>% 
            filter(pval.adj <= pval.cutoff) %>%
            mutate(source = src) %>%
            select(source, goid, ratioInSample, ratioInPop, pval.adj, gotype, goname)
        to = rbind(to, to1)
    }
    to
    #}}}
}

go_enrich_genesets <- function(tgs, pval.cutoff = 0.05) {
    #{{{
    te = tibble()
    for (tag1 in unique(tge$tag)) {
        gids = tge %>% filter(tag == tag1) %>% pull(gid)
        te1 = go_enrich_gosets(gids, pval.cutoff = pval.cutoff) %>% mutate(tag = tag1) %>%
            select(tag, everything())
        te = rbind(te, te1)
    }
    te
    #}}}
}

#fisher.test(matrix(c(hitInSample, hitInPop-hitInSample, sampleSize-hitInSample, failInPop-sampleSize +hitInSample), 2, 2), alternative='two.sided')




