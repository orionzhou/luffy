require(GenomicRanges)

f_gtb = file.path(DIR_Misc2, "genefam/31_merged.gtb")
g01 = read.table(f_gtb, header=T, as.is=T, sep="\t", quote="")
g01$strand[g01$strand==1] = "+"
g01$strand[g01$strand==-1] = "-"
loc_genes = GRanges(seqnames=g01$chr, ranges=IRanges(g01$beg, g01$end), strand=g01$strand, id=g01$id)

f_qry = file.path(DIR_Misc3, "spfinder/Mtruncatula/11_hmmsearch/07.bes")
q01 = read.table(f_qry, header=T, as.is=T, sep="\t", quote="")

locStr2Range <- function(str, seqid='chrN', type='') {
  tmp = gregexpr("complement([[:graph:]]+)", str)[[1]]
  strand = "+"
  if(tmp[1] == 1) { strand = "-" }
  tmp = gregexpr("([0-9]+\\.\\.[0-9]+)", str)[[1]]
  if(tmp[1] == -1) {
    GRanges()
  } else {
    strs = substring(str, tmp, tmp+attr(tmp, "match.length")-1)
    tmp = strsplit(strs, "\\.\\.")
    begs = as.numeric(sapply(tmp, "[", 1))
    ends = as.numeric(sapply(tmp, "[", 2))
    GRanges(seqnames=seqid, ranges=IRanges(begs, ends), strand=strand, type=rep(type, length(begs)))
  }
}
get_gene_model <- function(x) {
  seqid = as.character(x['chr'])
  beg = as.numeric(x['beg'])
  end = as.numeric(x['end'])
  strand = as.character(x['strand'])
  gr = GRanges(seqnames=seqid, ranges=IRanges(beg, end), strand=strand)
  elementMetadata(gr)[,"type"] = "all"

  grC = locStr2Range(x['loc_cds'], seqid, 'cds')
  grE = locStr2Range(x['loc_exon'], seqid, 'exon')
  
  grI = setdiff(gr, grE)
  if(length(grI) > 0) { elementMetadata(grI)[,"type"] = "intron" }

  gr5 = locStr2Range(x['loc_utr5'], seqid, 'utr5')
  gr3 = locStr2Range(x['loc_utr3'], seqid, 'utr3')
  
  GRangesList('all'=gr, 'cds'=grC, 'intron'=grI, 'utr5'=gr5, 'utr3'=gr3)
}
qry_gene_model <- function(gene, qLoc) {
  grC = intersect(qLoc, gene$cds)
  grI = intersect(qLoc, gene$intron)
  gr5 = intersect(qLoc, gene$utr5)
  gr3 = intersect(qLoc, gene$utr3)
  grO = setdiff(qLoc, gene$all)
  c('cds'=sum(width(grC)), 'intron'=sum(width(grI)), 'utr5'=sum(width(gr5)), 'utr3'=sum(width(gr3)), 'out'=sum(width(grO)))
}
qry_gene_models <- function(qry, genes_all) {
  seqid = as.character(qry['hId'])
  qLoc = locStr2Range(qry["hLoc"], seqid, qry["qId"])
  genes = subset(genes_all, chr==seqid & beg<=end(qLoc) & end>=start(qLoc))
  if(nrow(genes) == 0) {
    v = rep(NA, 6)
  } else {
    rownames(genes) = genes$id
    gLocs = apply(genes, 1, get_gene_model)
    df = data.frame(t(sapply(gLocs, qry_gene_model, qLoc=qLoc)))
    df = df[order(df$cds, decreasing=TRUE),]
    v = c(rownames(df)[1], as.numeric(df[1,]))
  }
  names(v) = c("id", "cds", "intron", "utr5", "utr3", "out")
  v
}

qrys = q01[1:10,]
r = apply(qrys, 1, qry_gene_models, genes_all = g01)
df = cbind(qrys[,c(-2,-7)], data.frame(t(r)))

chr = "chr5"
beg = 25793700
end = 25797499 
g11 = g01[g01$chr == chr & g01$beg <= end & g01$end >= beg,]
df = do.call("rbind", apply(g11, 1, get_loc_gene))

p <- ggplot(data=dummy, mapping=aes(x=x,y=y)) +
  facet_grid(par ~ ., scale='free') +
  layer(data=df1, geom='rect', mapping=aes(x=beg, y=0.5, xmin=beg, xmax=end, ymin=0, ymax=1, fill=type), geom_params=list(size=0)) +
  layer(data=df2, geom='rect', mapping=aes(x=beg, y=0.5, xmin=beg, xmax=end, ymin=1, ymax=2, fill=type), geom_params=list(size=0)) +
  scale_fill_manual(values=c('cds'='red', 'intron'='grey', 'utr5'='violet', 'utr3'='tan'), legend=TRUE) +
  labs(fill="Type") +
  scale_x_continuous(name='chr position (bp)', limits=c(beg, end), formatter='comma') +
  scale_y_continuous(name='', limits=c(0,10)) +
  theme_bw() +
  opts(axis.text.y=theme_blank())
ggsave(p, filename = file.path(dirW, "02_gd.png"), width=10, height=8)





