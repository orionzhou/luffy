f0 = file.path(DIR_Misc1, "window_stats/11_mt_35/03_gene.tbl")
gs = read.table(f0, header=TRUE, sep="\t")

tmp = apply(gs[,c(4:8)], 2, sum)
tmp2 = cbind(cat="genome", region=names(tmp), size=as.numeric(tmp))
tmp = apply(gs[gs$chr=='chr5',c(4:8)], 2, sum)
tmp3 = cbind(cat="chr5", region=names(tmp), size=as.numeric(tmp))
tmp = apply(d3[,c(4:8)], 2, sum)
tmp4 = cbind(cat="chr5_deletion", region=names(tmp), size=as.numeric(tmp))

d4 = rbind(tmp2,tmp3,tmp4)
d4 = data.frame(d4, stringsAsFactors=FALSE)

d4$cat = factor(d4$cat, levels=c("chr5_deletion", "chr5", "genome"))
d4$size = as.numeric(d4$size)
d4$region = factor(d4$region, levels=c("bp_cds", "bp_intron", "bp_utr5", "bp_utr3", "bp_intergenic"))

p = ggplot(d4) +
  geom_bar(aes(x=cat, y=size, fill=region), stat='identity', position='fill') +
  scale_fill_brewer(breaks=levels(d4$region), labels = c("CDS", "Intron", "UTR5", "UTR3", "Intergenic"), palette='Set1') +
  scale_y_continuous(name="", formatter='percent') +
  scale_x_discrete(name="") +
    opts(axis.text.x=theme_text(size=8), axis.text.y=theme_text(hjust=1, size=10, colour="blue")) +
  coord_flip()
ggsave(p, filename = file.path(dirW, "deletion04.png"), width=5, height=2);

a1 = readAssembly("mt_35")
a = drop.levels(subset(a1, type!='centromere'))
c = drop.levels(subset(a1, type=='centromere'))

f5 = file.path(DIR_Misc1, "window_stats/11_mt_35/05_gene_te.tbl")
d11 = read.table(f5, header=TRUE, sep="\t")
f6 = file.path(DIR_Misc1, "window_stats/11_mt_35/22_deletion.tbl")
d12 = read.table(f6, header=TRUE, sep="\t")
d13 = merge(d11, d12, by.x=c('chr','start','end'), by.y=c('chr','start','end'))

d14 = cbind(d13, length=d13$end-d13$start+1)
d15 = cbind(d14, percent_gene=d14$gene/d14$length, percent_TE=d14$transposable_element_gene/d14$length)
p = ggplot(data = d15) + 
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=percent_gene, fill='gene'), size=0) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=percent_gene, ymax=percent_gene+percent_TE, fill='TE'), size=0) +
  layer(data=a, geom='rect', mapping=aes(xmin=start, xmax=end, ymin=0.95, ymax=1, fill=type), geom_params=list(size=0)) +
  layer(data=c, geom='point', mapping=aes(x=(start+end)/2, y=0.9, shape=type), geom_params=list(size=2)) +
  scale_fill_manual(values=c('phase 1'='aquamarine', 'phase 3'='aquamarine4', 'gene'='violet', 'TE'='tan'), legend=TRUE) +
  scale_shape_manual(values=c('centromere'=17), legend=TRUE) +
  scale_x_continuous(name='chr position (bp)', formatter='comma') +
  scale_y_continuous(name='', limits=c(0,1), formatter="percent") +
  opts(axis.text.y=theme_text(hjust=1, size=8, colour="blue")) +
  facet_grid(chr ~ .) 
ggsave(p, filename = file.path(dirW, "deletion_dist.png"), width=10, height=8);

p = ggplot(data = d15) + 
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=cnt_deletion, fill='deletion'), size=0) +
  scale_shape_manual(values=c('centromere'=17), legend=TRUE) +
  scale_x_continuous(name='chr position (bp)', formatter='comma') +
  scale_y_continuous(name='number of deletion in 100kb windows') + 
  scale_fill_manual(values=c('deletion'='royalblue'), legend=TRUE) +
  opts(axis.text.y=theme_text(hjust=1, size=8, colour="blue")) +
  facet_grid(chr ~ .) 
ggsave(p, filename = file.path(dirW, "deletion_dist2.png"), width=10, height=8);

#genotyping
accs = get_mt_ids("acc26")
dir = file.path(DIR_Repo, "mt_35/40_sv/41_shared")
s1 = read.table(file.path(dir, "11.tbl"), sep="\t", as.is=T, header=T, quote="", comment.char="/")
gt = matrix(NA, nrow=nrow(s1), ncol=length(accs)+1)
colnames(gt) = c("id", accs)
for (i in 1:nrow(s1)) {
  id = s1$id_pindel[i]
  f_geno = file.path(dir, "66_exonerate", paste(id, ".tbl", sep=""))
  g1 = read.table(f_geno, sep="\t", as.is=T, header=T, quote="", comment.char="/")
  t1 = table( g1[,c('acc','support')] )
  gt[i, 'id'] = id
  for (acc in accs) {
    allele = "N"
    if(acc %in% rownames(t1)) {
      n_reads_ref = t1[acc, 'ref'] 
      n_reads_alt = t1[acc, 'alt'] 
      if(n_reads_ref > 2 && n_reads_alt <= 2) allele = "A"
      if(n_reads_ref <= 2 && n_reads_alt > 2) allele = "G"
    }
    gt[i, acc] = allele
  }
}
write.table(gt, file.path(dir, "71_genotype.tbl"), sep="\t", quote=F, row.names=F, col.names=T)
#manually transform to simpleSNP format
#sspFilter -i 01.ssp -o 02.ssp
#sspConvert -i 02.ssp -o 03.phy -f phylip

#seqret 03.phy 04.aln -osformat aln
#clustalw2 -INFILE=04.aln -BOOTSTRAP=1000 -OUTORDER=INPUT -OUTPUTTREE=phylip -BOOTLABELS=node -CLUSTERING=NJ 

#PhyML -i 03.phy -d nt

#calculate SFS
gt = read.table(file.path(dir, "71_genotype.tbl"), sep="\t", header=T, as.is=T)
cnt_allele <- function(x, allele) length(grep(allele, x))
n_N = apply(gt, 1, cnt_allele, allele="N")
n_ref = apply(gt, 1, cnt_allele, allele="A")
n_alt = apply(gt, 1, cnt_allele, allele="G")
sfs = data.frame(id=gt$id, n_N=n_N, n_ref=n_ref, n_alt=n_alt)
sfs.1 = sfs[sfs$n_ref>0 & sfs$n_alt>0, ]
sfs.2 = cbind(sfs.1, freq_alt=sfs.1$n_alt/(sfs.1$n_alt+sfs.1$n_ref))
write.table(sfs.2, file.path(dir, "74_sfs.tbl"), sep="\t", quote=F, row.names=F, col.names=T)

tmp = cut(sfs.2$freq_alt, breaks=c(seq(0,1,0.1)))
data = data.frame(DAF_bin=tmp)
p = ggplot(data) +
  geom_bar(aes(x=factor(DAF_bin)), width=0.4) + 
  scale_x_discrete(name="Proportion of Accessions with the Deletion") +
  scale_y_continuous(name="Number of Events") +
  opts(axis.text.x = theme_text(angle=45, size=8))
ggsave(file.path(dir, "74_sfs.png"), p, width=4, height=4)


#sv_effect.pl
dir = file.path(DIR_Repo, "mt_35/40_sv/61_effect")
s1 = read.table(file.path(dir, "01.tbl"), sep="\t", header=T, as.is=T)
e1 = read.table(file.path(dir, "02_gene.tbl"), sep="\t", as.is=T, header=T, quote="")
m1 = cbind(s1, e1[,-c(1:3)])
write.table(m1, file.path(dir, "03_merged.tbl"), sep="\t", quote=F, row.names=F, col.names=T)

ta1 = read.table(file.path(DIR_Misc2, "genefam", "41_genefam.tbl"), sep="\t", header=T, as.is=T, quote="")
ids_nbs = ta1$id[ta1$type1=='nbs_gene']
ids_sv_o = ids_sv_o = m1$id[m1$bp_cds>0]
ids_sv = unlist(strsplit(ids_sv_o, " "))
ids_sv[ids_sv %in% ids_nbs]


#plot deletion size distribution & SFS
dir = file.path(DIR_Repo, "mt_35/40_sv/90_stat")
m1 = read.table(file.path(dir, "../61_effect/03_merged.tbl"), sep="\t", as.is=T, header=T, quote="")
m2 = m1[,c('id_pindel', 'beg', 'end', 'size_d', 'size_i', 'bp_cds', 'bp_intron', 'bp_utr5', 'bp_utr3', 'bp_intergenic')]
s1 = read.table(file.path(dir, "../41_shared/74_sfs.tbl"), sep="\t", as.is=T, header=T, quote="")
data = merge(m2, s1, by.x='id_pindel', by.y='pos')

tmp = cut(data$size_d/1000, breaks=c(seq(0,5,0.2)))
p = ggplot(data.frame(size=tmp)) +
  geom_bar(aes(x=factor(size)), width=0.8) + 
  scale_x_discrete(name="Deletion Size (kb)") + 
  scale_y_continuous(name="Number Events") +
  opts(axis.text.x = theme_text(angle=45, size=8))
ggsave(file.path(dir, "01_size.png"), p, width=5, height=4)

tmp1 = cut(data$freq_der[data$bp_cds>0], breaks=c(seq(0,0.5,0.05)))
tmp2 = cut(data$freq_der[data$bp_cds==0], breaks=c(seq(0,0.5,0.05)))
tmp = rbind(data.frame(DAF_bin=tmp1, type='CDS'), data.frame(DAF_bin=tmp2, type='non-CDS'))
p = ggplot(tmp) +
  geom_bar(aes(x=factor(DAF_bin), fill=type), width=0.6, position='dodge') + 
  scale_fill_brewer(palette='Set1') +
  scale_x_discrete(name="Derived Alelle Frequency (folded)") +
  scale_y_continuous(name="Number of Deletion Events") +
  opts(axis.text.x = theme_text(angle=45, size=8))
ggsave(file.path(dir, "02_sfs.png"), p, width=5, height=4)


