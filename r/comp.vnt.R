require(rtracklayer)
require(GenomicRanges)
#require(PopGenome)
#require(pegas)
#require(VariantAnnotation)
require(ggplot2)
require(grid)
source("comp.fun.R")
source("Location.R")

dirw = file.path(Sys.getenv("misc3"), "comp.vnt")
diro = file.path(Sys.getenv("misc3"), "comp.stat")

##### obtain regions covered by 9 (out of 12) accessions
qnames = get_orgs(opt = "ingroup")
tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

gr = GenomicRanges::setdiff(grt, grp)
gra = GRanges()
for (qname in qnames) {
  diri = sprintf("%s/%s_%s", Sys.getenv("misc3"), qname, tname)
  
  fi = file.path(diri, '23_blat/31.9/gax')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:3,10)]
  colnames(ti) = c('chr', 'beg', 'end', 'lev')
  ti = ti[ti$lev <= 1,]
  grg = with(ti, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
  gra = c(gra, grg)
}
gr_cvg = as(coverage(gra), "GRanges")
for (i in 1:length(qnames)) {
  cat(sprintf("%2d: %d\n", i, sum(width(gr_cvg[mcols(gr_cvg)$score >= i, ]))))
}
grc = reduce(gr_cvg[mcols(gr_cvg)$score >= 9, ])

tcvg = data.frame(chr = seqnames(grc), beg = start(grc), end = end(grc))
fo = file.path(dirw, "81.cvg.tbl")
write.table(tcvg, fo, sep = "\t", row.names = F, col.names = F, quote = F)

##### obtain reference-mapping based covered regions (9 out of 12)
dirw = file.path(Sys.getenv("misc3"), "hapmap/12_ncgr")

qnames = get_orgs(opt = "ingroup")

gra = GRanges()
for (qname in qnames) {
  fi = sprintf("%s/35_cov/%s.bed", dirw, qname)
  ti = read.table(fi, header = F, sep = "\t", as.is = T)
  colnames(ti) = c('chr', 'beg', 'end')
  grg = with(ti, GRanges(seqnames = chr, ranges = IRanges(beg+1, end = end)))
  gra = c(gra, grg)
}
gr_cvg = as(coverage(gra), "GRanges")
for (i in 1:length(qnames)) {
  cat(sprintf("%2d: %d\n", i, sum(width(gr_cvg[mcols(gr_cvg)$score >= i, ]))))
}
grc = reduce(gr_cvg[mcols(gr_cvg)$score >= 9, ])

tcvg = data.frame(chr = seqnames(grc), beg = start(grc), end = end(grc))
fo = file.path(dirw, "81.cvg.tbl")
write.table(tcvg, fo, sep = "\t", row.names = F, col.names = F, quote = F)

##### misc tests of packages
#popgenome
dv = file.path(dir, "xxx")
x1 = readData(dv, include.unknown = T, format = "VCF")
get.sum.data(x1)

x1 = detail.stats(x1)
mafs = x1@region.stats@minor.allele.freqs[[1]]

x1 = F_ST.stats(x1)
x1 = diversity.stats(x1)
x1@nuc.diversity.within

# variantannotation + snpStats
require(VariantAnnotation)
require(snpStats)

fv = file.path(dir, "xxx/x.vcf")
vcf <- readVcf(fv, "mt40")
res = genotypeToSnpMatrix(vcf)
cat(nrow(res$map), "sites", sum(res$map$ignore), "poly-allelic\n", sep = " ")
snpsum = col.summary(res$genotypes)[!res$map$ignore,]
sum( (2*snpsum$Calls)/(2*snpsum$Calls-1) * 2 * snpsum$MAF * (1-snpsum$MAF) )

##### variant AFS
dirw = file.path(Sys.getenv("misc3"), "comp.vnt")
fa = file.path(dirw, "81.cvg.tbl")
ta = read.table(fa, header = F, sep = "\t", as.is = T)
gra = with(ta, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

fv = file.path(dirw, '25.stat.tbl')
tv = read.table(fv, header = T, sep = "\t", as.is = T)
gsnp = with(tv, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))

effmap = c(
  'downstream_gene_variant' = 'intergenic',
  'upstream_gene_variant' = 'intergenic',
  'intergenic_region' = 'intergenic',
  'missense_variant' = 'replacement',
  'synonymous_variant' = 'synonymous',
  'stop_retained_variant' = 'synonymous',
  'initiator_codon_variant' = 'synonymous',
  'intron_variant' = 'intron',
  'splice_region_variant' = 'intron',
  '5_prime_UTR_variant' = 'utr',
  '3_prime_UTR_variant' = 'utr',
  '5_prime_UTR_premature_start_codon_gain_variant' = 'utr',
  'stop_gained' = 'large_effect',
  'splice_acceptor_variant' = 'large_effect',
  'splice_donor_variant' = 'large_effect',
  'stop_lost' = 'large_effect',
  'start_lost' = 'large_effect'
)
ff = file.path(dirw, "23.snp.anno.tbl")
tf = read.table(ff, header = F, sep = "\t", as.is = T)
stopifnot(nrow(tv) == nrow(tf))
stopifnot(sum(!tf$V3 %in% names(effmap)) == 0)


to = cbind(tv, eff = as.character(effmap[tf$V3]))
to = to[to$nsam >= 5 & to$nalt == 1,]
to = cbind(to, altcnt = round(to$nsam * to$aaf))
to = cbind(to, mac = pmin(to$altcnt, to$nsam-to$altcnt))
to = cbind(to, maf = to$mac / to$nsam)

intvs = c(seq(0, 0.42, by = 0.07), 0.5)
labs = intvs[-1]
to = cbind(to, maf.bin = cut(to$maf, intvs, include.lowest = T))
#to = cbind(to, aaf.bin = cut(to$aaf, intvs, include.lowest = T))

do1 = ddply(to, .(eff), summarise, eff_tot_cnt = length(eff))
do2 = ddply(to, .(eff, maf.bin), summarise, eff_cnt = length(eff))
do = merge(do2, do1, by = 'eff')
do = cbind(do, prop = do$eff_cnt / do$eff_tot_cnt)

fo = sprintf("%s/16_snp_maf.pdf", diro)
p1 = ggplot(do) +
  geom_bar(mapping = aes(x = maf.bin, y = prop, fill = eff), 
    stat = 'identity', position = 'dodge', 
    geom_params = list(width = 0.8, alpha = 0.8)) +
  scale_fill_brewer(palette='Dark2', name = '') +
  scale_x_discrete(name = 'Minor Allele Frequency') +
  scale_y_continuous(name = 'Proportion', expand = c(0,0), limits = c(0,0.52)) +
  theme_bw() +
  theme(legend.position = c(0.7, 0.8), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.7, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_blank(), legend.text = element_text(size = 8, angle = 0)) +
  theme(plot.margin = unit(c(0,1,1,0), "lines")) +
  theme(axis.text.x = element_text(size = 8, colour = "grey", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0))
ggsave(p, filename = fo, width = 5, height = 5)


fv = sprintf("%s/52.stat.tbl", dirw)
tv = read.table(fv, sep = '\t', header = T, as.is = T)
tv = cbind(tv, size = (tv$rsize + tv$asize - 2), altcnt = round(tv$nsam * tv$aaf))
tv$nsam = tv$nsam + 1
tv = tv[tv$nsam >= 5 & tv$nalt == 1,]
tv = cbind(tv, type = NA)
tv$type[tv$size < 50 & tv$asize == 1] = 'small del'
tv$type[tv$size < 50 & tv$rsize == 1] = 'small ins'
tv$type[tv$size < 50 & tv$rsize > 1 & tv$asize > 1] = 'small mnp'
tv$type[tv$size >= 50 & tv$asize == 1] = 'large del'
tv$type[tv$size >= 50 & tv$rsize == 1] = 'large ins'
tv$type[tv$size >= 50 & tv$rsize > 1 & tv$asize > 1] = 'large mnp'

to = tv
colnames(to)[ncol(to)] = 'eff'
to = to[to$nsam >= 5 & to$nalt == 1,]
to = cbind(to, mac = pmin(to$altcnt, to$nsam-to$altcnt))
to = cbind(to, maf = to$mac / to$nsam)

intvs = c(seq(0, 0.42, by = 0.07), 0.5)
labs = intvs[-1]
to = cbind(to, maf.bin = cut(to$maf, intvs, include.lowest = T))
#to = cbind(to, aaf.bin = cut(to$aaf, intvs, include.lowest = T))

do1 = ddply(to, .(eff), summarise, eff_tot_cnt = length(eff))
do2 = ddply(to, .(eff, maf.bin), summarise, eff_cnt = length(eff))
do = merge(do2, do1, by = 'eff')
do = cbind(do, prop = do$eff_cnt / do$eff_tot_cnt)

fo = sprintf("%s/16_sv_maf.pdf", diro)
p2 = ggplot(do) +
  geom_bar(mapping = aes(x = maf.bin, y = prop, fill = eff), 
    stat = 'identity', position = 'dodge', 
    geom_params = list(width = 0.8, alpha = 0.8)) +
  scale_fill_brewer(palette='Dark2', name = '') +
  scale_x_discrete(name = 'Minor Allele Frequency') +
  scale_y_continuous(name = 'Proportion', expand = c(0,0), limits = c(0,0.8)) +
  theme_bw() +
  theme(legend.position = c(0.7, 0.8), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.7, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_blank(), legend.text = element_text(size = 8, angle = 0)) +
  theme(plot.margin = unit(c(0,1,1,0), "lines")) +
  theme(axis.text.x = element_text(size = 8, colour = "grey", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0))
ggsave(p2, filename = fo, width = 5, height = 5)


##### SNP ThetaPi ThetaW for cds/intron/intergenic/utr/syn/nonsyn
fa = file.path(dirw, "81.cvg.tbl")
ta = read.table(fa, header = F, sep = "\t", as.is = T)
gra = with(ta, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

f1 = file.path(Sys.getenv("genome"), "HM101", "51.tbl")
t1 = read.table(f1, header = F, sep = "\t", as.is = T)
colnames(t1) = c("chr", "beg", "end", "srd", "id", "type", "fam")

tt = t1[t1$type == 'cds',]
gcds = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gcds = GenomicRanges::intersect(gra, reduce(gcds))

tt = t1[t1$type == 'intron',]
gito = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gito = GenomicRanges::intersect(GenomicRanges::setdiff(gra, gcds), reduce(gito))

tt = t1[t1$type == 'utr5',]
gut5 = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gut5 = GenomicRanges::intersect(GenomicRanges::setdiff(gra, c(gcds, gito)), reduce(gut5))

tt = t1[t1$type == 'utr3',]
gut3 = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gut3 = GenomicRanges::intersect(GenomicRanges::setdiff(gra, c(gcds, gito, gut5)), reduce(gut3))

gitr = GenomicRanges::setdiff(gra, reduce(c(gcds, gito, gut5, gut3)))

# reading the huge file is very slow
f2 = file.path(Sys.getenv("genome"), "HM101", "51.codon.degen.tbl")
t2 = read.table(f2, header = F, sep = "\t", as.is = T)
gsyn = with(t2[t2$V3 == 4,], GRanges(seqnames = V1, ranges = IRanges(V2, end = V2)))
gsyn = GenomicRanges::intersect(gra, reduce(gsyn))
grpl = with(t2[t2$V3 == 0,], GRanges(seqnames = V1, ranges = IRanges(V2, end = V2)))
grpl = GenomicRanges::intersect(gra, reduce(grpl))

grl = GRangesList(cds = gcds, synonymous = gsyn, replacement = grpl, intron = gito, utr5 = gut5, utr3 = gut3, intergenic = gitr)

  pi = sum(tv$nucdiv) / sum(width(gra))
  tw = table(tv$nsam)
  ans = sapply(as.numeric(names(tw)), get_har <- function(x) sum(1 / 1:(x-1)))
  thetaw = sum(as.numeric(tw) / ans) / sum(width(gra))

stats = list()
stats[["Total"]] = matrix(c(
    format(sum(width(gra)), big.mark = ","), 
    '-', 
    format(nrow(tv), big.mark = ","),
    sprintf("%.04f", pi), 
    sprintf("%.04f", thetaw)), nrow = 1, 
    dimnames = list(NULL, c("Covered bases", "Percent",
    "Polymorphic sites", "Pi", "ThetaW")))
for (i in 1:length(grl)) {
  ma = as.matrix(findOverlaps(gsnp, grl[i]))
  idxs = ma[,1]
  tvs = tv[idxs,]
  
  bp = sum(width(grl[[i]]))
  pct = bp / sum(width(gra))
  n_snp = nrow(tvs)
  pi = sum(tvs$nucdiv) / bp
  
  tw = table(tvs$nsam)
  ans = sapply(as.numeric(names(tw)), get_har <- function(x) sum(1 / 1:(x-1)))
  thetaw = sum(as.numeric(tw) / ans) / bp

  stat = c(
    format(bp, big.mark = ","), 
    sprintf("%.02f", pct), 
    format(n_snp, big.mark = ","),
    sprintf("%.04f", pi), 
    sprintf("%.04f", thetaw))
  stats[[names(grl)[i]]] = matrix(stat, nrow = 1, 
    dimnames = list(NULL, c("Covered bases", "Percent",
    "Polymorphic sites", "Pi", "ThetaW")))
}
ds = do.call(rbind.data.frame, stats)
fo = file.path(diro, "13_diversity.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)

