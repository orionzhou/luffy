require(GenomicRanges)

f_gs = file.path(DIR_Data, "genome/mt_35/51_size.tbl")
gs = read.table(f_gs, sep="\t", as.is=T, header=T)

f_gff = file.path(DIR_Data, "genome/mt_35/10_model/55.tbl")
g01 = read.table(f_gff, sep="\t", as.is=T, header=T, quote="")

#g01.bak = g01
#g01 = g01.bak[490006:500005,]

g02 = g01[g01$type == 'gene' | g01$type == 'transposable_element_gene',]
gene = GRanges(seqnames=g02$chr, ranges=IRanges(g02$beg, g02$end), strand=g02$strand, type=g02$type)
names(gene) = g02$Name
seqlengths(gene)[gs$chr] = gs$length_total

df2Grl <- function(df, type) {
  if(type == 'transcript') {
    gr = GRanges(seqnames=df$chr, ranges=IRanges(df$beg, df$end), strand=df$strand, conf_class=df$conf_class, Note=df$Note)
    names(gr) = df$Name
  } else if(type == 'exon') {
    gr = GRanges(seqnames=df$chr, ranges=IRanges(df$beg, df$end), strand=df$strand)
  } else if(type == 'cds') {
    gr = GRanges(seqnames=df$chr, ranges=IRanges(df$beg, df$end), strand=df$strand, phase=df$phase)
  } else {
    stop(paste("unsupported type: ", type, sep=""))
  }
  gr
}

g03 = g01[g01$type == 'mRNA',]
transcript = GRangesList(dlply(g03, .(Parent), df2Grl, type='transcript'))
seqlengths(transcript)[gs$chr] = gs$length_total

g04 = g01[g01$type == 'exon',]
exon = GRangesList(dlply(g04, .(Parent), df2Grl, type='exon'))
seqlengths(exon)[gs$chr] = gs$length_total

g05 = g01[g01$type == 'CDS',]
cds = GRangesList(dlply(g05, .(Parent), df2Grl, type='cds'))
seqlengths(cds)[gs$chr] = gs$length_total

save(gene, transcript, exon, cds, file=file.path(DIR_RData, "mt35_gene.RData"))
#load(file.path(DIR_RData, "mt35_gene.RData"))
