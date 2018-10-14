require(tidyverse)
require(Hmisc)

dirw = '/home/springer/zhoux379/data/misc1/maize.acr'

### read in bed file
fb = file.path(dirw, "01.bed")
tb = read_tsv(fb, col_names = T, col_types = 'ciiidc') %>%
    mutate(start = start + 1,
           size = stop - start + 1,
           lid = sprintf("%s-%d-%d",chr, start, stop))
tb
length(unique(tb$lid))
colnames(tb)[5:6] = c("gid", "gdist")
length(unique(tb$gid))

describe(tb$size)

table(tb$category)
#describe(tb$gdist[tb$category == 'proximal']) # 1-2000
#describe(tb$gdist[tb$category == 'distal']) # 2000-850000
#describe(tb$gdist[tb$category == 'intra']) # 0

### add PH207 synteny genes
fg = '/home/springer/zhoux379/data/genome/Zmays_v4/51.gtb'
tg = read.table(fg, sep = "\t", header = T, as.is = T)

fm = '/home/springer/zhoux379/data/genome/PH207/synteny/mapping.tsv'
tm = read_tsv(fm, col_names = F)
colnames(tm) = c("bchr","bbeg",'bend','bgid','pchr','pbeg','pend','pgid','type1','orthoType')
tm2 = tm  %>% select(bgid,pgid,pchr,pbeg,pend,orthoType)
idxs = which(tm2$pchr %in% sprintf("%d", 1:10))
tm2$pchr[idxs] = sprintf("chr%02d", as.integer(tm2$pchr[idxs]))
tm3 = unique(tm2)
nrow(tm3) - nrow(unique(tm3[,1:2]))

tb2 = tb %>% left_join(tm3, by = c('gid' = 'bgid'))

### align to PH207 and process outputs
# cd $misc1/jackie.acr
# seqret.pl -d $genome/Zmays_v4/11_genome.fas -b 01.bed -o 02.fas
# blat $genome/PH207/21.blat/db.2bit 02.fas -ooc=$genome/PH207/21.blat/db.2bit.tile11.ooc 04.psl
# psl2tsv.pl -i 04.psl -o 05.tsv

fi = file.path(dirw, "05.tsv")
ti = read_tsv(fi)[,1:20]
#ti = ti[ti$alnLen/ti$qSize >= 0.8 & ti$ident >= 0.8,]

tb4 = tb2 %>%
    left_join(ti, by = c('lid' = 'qId')) %>%
	mutate(pgdist = ifelse(is.na(pchr) | tId != pchr, NA,
		ifelse((tBeg >= pbeg & tBeg <= pend) | (tEnd >= pbeg & tEnd <= pend), 0, 
		ifelse(tEnd < pbeg, pbeg - tEnd - 1, 
		ifelse(tBeg > pend, tBeg - pend - 1, NA)))))
tb4$pgdist[!is.na(tb4$pgdist) & tb4$pgdist > 1000000] = NA
table(is.na(tb4$pgdist))

tb41 = tb4[!is.na(tb4$pgdist),]
tb42 = tb4[is.na(tb4$pgdist),]
length(unique(tb41$lid))

# ACRs with syntenic hits
tb41s = tb41 %>% 
    group_by(lid, pchr, pbeg, pend) %>%
    arrange(-score, pgdist) %>%
    summarise(nbest = sum(score == score[1] & pgdist == pgdist[1]))
table(tb41s$nbest)

to = tb41 %>% 
    group_by(lid, pchr, pbeg, pend) %>%
    #group_by(lid, gid, chr, start, stop, summit, gdist, category, size, pgid, pchr, pbeg, pend, orthoType) %>%
    arrange(-score, pgdist) %>%
    slice(1) %>%
    mutate(varType = 
    ifelse(alnLen == qSize & misMatch == 0 & qNumIns+tNumIns == 0, 'Identical', 
    ifelse(alnLen < qSize & misMatch == 0 & qNumIns+tNumIns == 0, 'assemblyGap',
    ifelse(misMatch > 0 & qNumIns+tNumIns == 0, 'onlySNP', 
    ifelse(qBaseIns >= 100 & qNumIns == 1, "qIns100",
    ifelse(tBaseIns >= 100 & tNumIns == 1, "tIns100", "Indel"))
    ))))
table(to$varType)



outputs = c(
'',
sprintf("%5d total ACRs:", length(unique(tb4$lid))),
sprintf("%5d do not map to PH207", length(unique(tb4$lid[is.na(tb4$tId)]))),
sprintf("%5d mapped to PH207:", length(unique(tb4$lid[!is.na(tb4$tId)]))),
sprintf("   %5d mapped in synteny (1Mb within the orthologous PH207 gene)", length(unique(tb4$lid[!is.na(tb4$pgdist)]))),
sprintf("       %5d full conservation (100%% coverage, no SNP, no InDel)", sum(to$varType=='Identical')),
sprintf("       %5d likely near assembly gaps (100%% coverage, no SNP, no InDel)", sum(to$varType=='assemblyGap')),
sprintf("       %5d only SNP(s), no Indel", sum(to$varType=='onlySNP')),
sprintf("       %5d InDel(s) and/or SNPs", sum(to$varType=='Indel')),
sprintf("       %5d contains a B73 insertion of >=100bp", sum(to$varType=='qIns100')),
sprintf("       %5d contains a PH207 insertion of >=100bp", sum(to$varType=='tIns100')),
''
)
cat(paste(outputs, collapse = "\n"))


#to[,c('chr','start','stop','summit','gid','gdist','category', 
    #'qId','qBeg','qEnd','qSrd','qSize',
    #'tId','tBeg','tEnd','tSrd','tSize',
    #'alnLen','match','misMatch','baseN',
    #'qNumIns','tNumIns','qBaseIns','tBaseIns','ident',
    #'pgid','pgdist','varType')]
tp = to %>% transmute(chr = chr, start = start, stop = stop, summit = summit, 
                      gid = gid, gdist = gdist, category = category,
                      match = match, misMatch = misMatch, 
                      numIns = qNumIns+tNumIns, baseIns = qBaseIns + tBaseIns,
                      pgid = pgid, pgdist = pgdist, varType = varType)
fo = file.path(dirw, "09.tsv")
write.table(tp, fo, sep = "\t", row.names = F, col.names = T, quote = F)

for (type in c("Identical", "qIns100", "tIns100")) {
    fo = sprintf("%s/09.%s.tsv", dirw, type)
    write.table(tp[tp$varType == type,], fo, sep = "\t", row.names = F, col.names = T, quote = F)
}
