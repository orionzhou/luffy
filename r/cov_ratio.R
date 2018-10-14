plotCov <- function() {
  db = "mt_30";
  a = readAssembly(db);
  f = file.path(DIR_Misc2, 'cov', '11acc_10kb.txt');
  cov <- read.table(f, sep = "\t", header = TRUE);

  acc = "HM005";
  cov2 <- cov[cov$acc==acc, ];
  p <- ggplot(data=cov2) + 
    geom_rect(mapping=aes(xmin=start, xmax=start+10000, ymin=cov1, ymax=cov2), size=0, alpha=0.8) +
    facet_grid(chr ~ .) + 
    scale_colour_brewer(palette="Set1") +
    scale_x_continuous(name='chr position (bp)', formatter='comma') +
    scale_y_continuous(trans='log2') +
    opts(legend.title=theme_blank());
  ggsave(p, filename = file.path(DIR_Stat, paste(acc, ".png", sep="")), width=15, height=8);

  accs = c("HM008", "HM006");
  cov.1 = cov[cov$acc==accs[1],];
  cov.2 = cov[cov$acc==accs[2],];
  covR = cbind(cov.1[,c(1:3)], cov.1[,c(6:7)], cov.2[,c(6:7)], NA, NA);
  colnames(covR)[4:9]=c('cov1', 'uniqcov1', 'cov1', 'uniqcov2', 'covr', 'uniqcovr');
  ratioCorrection = sum(covR$uniqcov1) / sum(covR$uniqcov2);
  covR$uniqcovr = log2((covR$uniqcov1 / covR$uniqcov2) / ratioCorrection);
  p <- ggplot(data=covR) + 
    geom_point(mapping=aes(x=start, y=uniqcovr), geom_params=list(size=0.7, alpha=0.8)) +
#    layer(data=l, geom='point', mapping=aes(x=pos, y=0, color=flag), position=position_jitter(w=2, h=4), geom_params=list(size=0.8, alpha=1)) +
    layer(data=a, geom='rect', mapping=aes(xmin=start, xmax=end, ymin=5.5, ymax=6, fill=type), geom_params=list(size=0)) +
    scale_fill_manual(values=c('phase 1 BAC'='orchid1', 'phase 2 BAC'='salmon', 'phase 3 BAC'='red1', 'centromere'='black', 'BAC of unknown phase'='green'), legend=TRUE) +
    scale_colour_brewer(palette="Set1") +
    scale_x_continuous(name='chr position (bp)', formatter='comma') +
    scale_y_continuous(name='log2(ratio)', limits=c(-6,6)) +
    facet_grid(chr ~ .) + 
    opts(legend.title=theme_blank(), title=paste('unique coverage ratio[', accs[1], ":", accs[2], "]", sep=""));
  ggsave(p, filename = file.path(DIR_Stat, "cov_ratio", paste(accs[1], "-", accs[2], ".png", sep="")), width=15, height=8);
}

plot_cov_info <- function(cov_info) {
  #plot total coverage / unique coverage per accession
  p = ggplot(data = cov_info) + 
    geom_bar(mapping = aes(x=chr, y=cov1_sum/pos_sum, fill=chr), stat='identity', position='dodge') +
    scale_fill_brewer(palette='Set1') +
    facet_wrap( ~ acc, nrow = 6) +
    scale_x_discrete(name="") +
    scale_y_continuous(name="coverage") +
    opts(axis.text.x = theme_blank());
  ggsave(p, filename=file.path(DIR_Misc1, "r/cov", "01_cov.png"), width=12, height=12)
  p = ggplot(data = cov_info) + 
    geom_bar(mapping = aes(x=chr, y=cov2_sum/pos_sum, fill=chr), stat='identity', position='dodge') +
    scale_fill_brewer(palette='Set1') +
    facet_wrap( ~ acc, nrow = 6) +
    scale_x_discrete(name="") +
    scale_y_continuous(name="unique coverage") +
    opts(axis.text.x = theme_blank());
  ggsave(p, filename=file.path(DIR_Misc1, "r/cov", "02_cov_uniq.png"), width=12, height=12)
}

fi = file.path(DIR_Repo, "mt_35", "30_vnt/21_coverage/11_info.tbl")
cov_info = read.table(fi, header=TRUE, sep="\t")
cov_info = cbind(cov_info, cov1_mean_chr=cov_info$cov1_sum/cov_info$pos_sum, cov2_mean_chr=cov_info$cov2_sum/cov_info$pos_sum)

acc1 = "HM018" #"HM005"
acc2 = "HM101" #"HM014"
cov_sum1 = cov_info[cov_info$acc==acc1,]
cov_sum2 = cov_info[cov_info$acc==acc2,]
ratio_adjust = data.frame(chr=cov_sum1$chr, adjust_cov1=cov_sum1$cov1_sum / cov_sum2$cov1_sum, adjust_cov2=cov_sum1$cov2_sum / cov_sum2$cov2_sum)

#coverage for aCGH probes
l1 = read.table(file.path(DIR_Misc1, "seq11/02_loc.tbl"), header=TRUE, sep="\t")
v1 = read.table(file.path(DIR_Misc1, "seq11/35_vnt_sum.tbl"), header=TRUE, sep="\t")
c1 = read.table(file.path(DIR_Misc1, "seq11/41_cov.tbl"), header=TRUE, sep="\t")

d1 = merge(v1, c1, by=c("id", "acc"), all.x=TRUE, all.y=TRUE)
d1$vnt_cnt[is.na(d1$vnt_cnt)] = 0
d2= merge(d1, l1, by.x='id', by.y='id', all.x=FALSE, all.y=FALSE)
cov = d2[, c(-3,-5,-6)]
colnames(cov)[12] = 'chr'

write.table(cov[,c(1, 12:14, 2:11)], file.path(DIR_R, "acgh/21_cov.tbl"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

cov_acc1 = cov[cov$acc==acc1,]
cov_acc2 = cov[cov$acc==acc2,]
cov_acc12 = merge(cov_acc1, cov_acc2, by=c('chr', 'beg', 'end'))
cov_acc12 = merge(cov_acc12, ratio_adjust, by='chr')
covr = data.frame(cov_acc12, cov2_ratio=log2((cov_acc12$cov2_mean.x / cov_acc12$cov2_mean.y) / cov_acc12$adjust_cov2))

covr2 = covr[,c(4,1:3,7,5:6,8:14,16:17,19:25,26:28)]
colnames(covr2)[c(1,5)] = c('id', 'length')
fname = paste(acc1, "_", acc2, ".tbl", sep='')
write.table(covr2, file.path(DIR_R, "acgh", fname), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

p = ggplot(data = covr2) + 
  geom_point(aes(x=(beg+end)/2000000, y=cov2_ratio), size=0.3) + 
  geom_hline(yintercept=c(0), color='purple', size=0.1) +
  geom_hline(yintercept=c(-1,1), color='blue', size=0.1) +
  geom_hline(yintercept=c(-2,2), color='red', size=0.1) +
  scale_x_continuous(name="chr position (Mbp)") +
  scale_y_continuous(name="log2(coverage ratio)", limits=c(-5,5)) +
  opts(axis.text.y=theme_text(size=8)) + 
  facet_grid(chr ~ .)
fname = paste(acc1, "_", acc2, ".png", sep='')
ggsave(p, filename=file.path(DIR_R, "acgh", fname), width=10, height=4)

#whole genome coverage ratio
fi = file.path(DIR_Misc2, "cov", "10k.tbl")
cov_10k = read.table(fi, header=TRUE, sep="\t")
fi = file.path(DIR_Misc2, "cov", "20k.tbl")
cov_20k = read.table(fi, header=TRUE, sep="\t")
fi = file.path(DIR_Misc2, "cov", "50k.tbl")
cov_50k = read.table(fi, header=TRUE, sep="\t")
fi = file.path(DIR_Misc2, "cov", "100k.tbl")
cov_100k = read.table(fi, header=TRUE, sep="\t")

winsize = "10k"
if(winsize=="10k") {cov=cov_10k} else if(winsize=="20k") {cov=cov_20k} else if(winsize=="50k") {cov=cov_50k} else if(winsize=="100k") {cov=cov_100k}
cov_acc1 = cov[cov$acc==acc1,]
cov_acc2 = cov[cov$acc==acc2,]
cov_acc12 = merge(cov_acc1, cov_acc2, by=c('chr', 'beg', 'end'))
cov_acc12 = merge(cov_acc12, ratio_adjust, by='chr')
covr = data.frame(cov_acc12[, -c(4,11)], cov2_ratio=log2((cov_acc12$cov2_mean.x / cov_acc12$cov2_mean.y) / cov_acc12$adjust_cov2))

dirW = file.path(DIR_R, "cov", paste(acc1, acc2, sep="_"))
plot_cov_lg(covr$cov2_mean.x, covr$cov2_mean.y, acc1, acc2, winsize, dirW)
plot_cov_dist(covr, winsize, dirW)

cutoff_cov2 = 2
for (i in 1:nrow(covr)) {
  if(covr$cov2_mean.x[i] < cutoff_cov2 || covr$cov2_mean.y[i] < cutoff_cov2) {
    covr$cov2_ratio[i] = NA
  }
}

plot_cov_dist <- function(covr, winsize, dirW) {
  p = ggplot(data = covr) + 
    geom_point(aes(x=(beg+end)/2000000, y=cov2_ratio), size=0.5) + 
    geom_hline(yintercept=c(0), color='purple', size=0.1) +
    geom_hline(yintercept=c(-1,1), color='blue', size=0.1) +
    geom_hline(yintercept=c(-2,2), color='red', size=0.1) +
    scale_x_continuous(name="chr position (Mbp)") +
    scale_y_continuous(name="log2(coverage ratio)", limits=c(-4,4)) +
    opts(axis.text.y=theme_text(size=8)) + 
    facet_grid(chr ~ .)
  ggsave(p, filename=file.path(dirW, paste(winsize, ".png", sep='')), width=10, height=5)
}

plot_cov_lg <- function(a, b, name1, name2, winsize, dirW) {
  fname = file.path(dirW, paste(winsize, "_lg.png", sep=''))
  png(fname)
  plot(a, b, type="p", xlab=name1, ylab=name2, xlim=c(0,80), ylim=c(0,80))
  title(main=winsize)
  fit = lm(b~a)
  abline(fit, col="blue")
  fit.sum = summary(fit)
  ann = paste("adjusted Rsquare = ", sprintf("%.04f", fit.sum$adj.r.squared), sep="")
  text(50, 0, ann, col='red')
  dev.off()
}

