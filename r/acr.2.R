require(tidyverse)
require(Hmisc)

dirw = '/home/springer/zhoux379/data/misc1/maize.acr'

# read in BED file
fb = file.path(dirw, "00.txt")
tb = read_tsv(fb, col_names = T, col_types = 'ciiidc') %>%
    mutate(start = as.integer(start + 1),
           size = stop - start + 1,
           lid = sprintf("%s-%d-%d",chr, start, stop))
tb
length(unique(tb$lid))
describe(tb$size)
table(tb$category)

### align to PH207 and process outputs - check acr.1.mapping.md

# read blat alignment, filter
fi = file.path(dirw, "05.tsv")
ti = read_tsv(fi)[,1:20] %>%
    filter(alnLen/qSize >= 0.8, ident >= 0.8)

# take best mappings and report alignment statistics
ti2 = ti %>%
    group_by(qId) %>%
    arrange(-score) %>%
    slice(1) %>%
    mutate(varType = 
    ifelse(alnLen == qSize & misMatch == 0 & qNumIns+tNumIns == 0, 'Identical', 
    ifelse(alnLen < qSize & misMatch == 0 & qNumIns+tNumIns == 0, 'assemblyGap',
    ifelse(misMatch > 0 & qNumIns+tNumIns == 0, 'onlySNP', "SNP+Indel"))))

ti3 = ti %>%
    group_by(qId) %>%
    arrange(-score) %>%
    summarise(bestMappings = sum(score==score[1]))

ti4 = ti2 %>% left_join(ti3, by = 'qId')

# join bed file with blat alignment
to = tb %>%
    left_join(ti4, by = c('lid' = 'qId')) %>%
    transmute(b.chrom = chr, b.start = start, b.stop = stop, b.size = size,
              summit = summit, count = count, category = category,
              bestMappings = bestMappings,
              p.chrom = tId, p.start = tBeg, p.stop = tEnd, varType = varType,
              alnLen = alnLen, match = match, misMatch = misMatch, 
              numIns = qNumIns+tNumIns, baseIns = qBaseIns+tBaseIns) %>%
    replace_na(list(bestMappings = 0, varType = 'unknown'))
table(to$varType, useNA = 'ifany')
table(to$bestMappings, useNA = 'ifany')

# print a text report
outputs = c(
'',
sprintf("%5d total ACRs:", nrow(to)),
sprintf("%5d did not map (to PH207)", sum(is.na(to$p.chrom))),
sprintf("%5d mapped (to PH207):", sum(!is.na(to$p.chrom))),
sprintf("    %5d full conservation (100%% coverage, no SNP, no InDel)", sum(to$varType=='Identical')),
sprintf("    %5d likely near assembly gaps (100%% coverage, no SNP, no InDel)", sum(to$varType=='assemblyGap')),
sprintf("    %5d only SNP(s), no Indel", sum(to$varType=='onlySNP')),
sprintf("    %5d InDel(s) and/or SNPs", sum(to$varType=='SNP+Indel')),
''
)
cat(paste(outputs, collapse = "\n"))

# write to tsv file
fo = file.path(dirw, "09.tsv")
write_tsv(to, fo)

