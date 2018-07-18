readIn <- function(fName) {
  d <- read.table(file.path(DIR_Stat, fName), header=TRUE);
  dl <-length(d);
  colnames(d)=c('geneFamily', 'geneId', 'type', "region", "acc", 'totCov', 
    'totVnt', 'calledVnt', 'coveredVnt', 'length');
  return (d);
}
agg1 <- function(d) {
  #aggregate by geneFamily & geneId & type & acc
  q=aggregate(d$totCov,by=list(d$geneFamily,d$geneId,d$type,d$acc),FUN=sum);
  w=aggregate(d$totVnt,by=list(d$geneFamily,d$geneId,d$type,d$acc),FUN=sum);
  e=aggregate(d$calledVnt,by=list(d$geneFamily,d$geneId,d$type,d$acc),FUN=sum);
  r=aggregate(d$coveredVnt,by=list(d$geneFamily,d$geneId,d$type,d$acc),FUN=sum);
  j=aggregate(d$length,by=list(d$geneFamily,d$geneId,d$type,d$acc),FUN=sum);  
  t=cbind(q,w[5],e[5],r[5],j[5]);
  tl <- length(t);
  t[tl+1]=t[5]/t[6];
  t[tl+1][is.na(t[tl+1])] <- 0;
  t[tl+2]=t[7]/(t[9]/t[5] * t[8]) * 10000;
  t[tl+2][is.na(t[tl+2])] <- 0;
  colnames(t)=c('geneFamily', 'geneId', 'type', "acc", 'totCov', 
    'totVnt', 'calledVnt', 'coveredVnt', 'length', 'avgCov', 'vntDensity');
  return (t);
}
plotT <- function(t) {
  geneNames = levels(t$geneFamily);
  genicRegion = levels(t$type);
  for (i in 1:length(geneNames)) {
    tt=t[t$geneFamily == geneNames[i],];
    title = geneNames[i];
    p <- ggplot(tt, aes(x=avgCov, fill=type)) + 
      xlim(c(0,70)) + xlab("Average Coverage") + 
      # geom_histogram(aes(y=..density..), binwidth=1, position="dodge") + ylim(c(0,50)) + 
      geom_density(alpha=0.2) + ylim(c(0,0.15)) + 
      facet_wrap(~ acc) +
      opts(title=title);
    ggsave(p, filename = file.path(DIR_Stat,paste("avgCov_",title,".png",sep="")), 
      width=8, height=5);
  }
}
agg2 <- function(d) {
  #aggregate by geneFamily & geneId & acc
  q=aggregate(d$totCov,by=list(d$geneFamily,d$geneId,d$acc),FUN=sum);
  w=aggregate(d$totVnt,by=list(d$geneFamily,d$geneId,d$acc),FUN=sum);
  e=aggregate(d$calledVnt,by=list(d$geneFamily,d$geneId,d$acc),FUN=sum);
  r=aggregate(d$coveredVnt,by=list(d$geneFamily,d$geneId,d$acc),FUN=sum);  
  j=aggregate(d$length,by=list(d$geneFamily,d$geneId,d$acc),FUN=sum);  
  g=cbind(q,w[4],e[4],r[4],j[4]);
  gl = length(g);
  g[gl+1]=g[4]/g[5];
  g[gl+1][is.na(g[gl+1])] <- 0;
  g[gl+2]=g[6]/(g[8]/g[5] *g[7]) * 10000;
  g[gl+2][is.na(g[gl+2])] <- 0;
  colnames(g)=c('geneFamily', 'geneId', "acc", 'totCov', 
    'totVnt', 'calledVnt', 'coveredVnt', 'length', 'avgCov', 'vntDensity');
}
plotG <- function(g) {
  p <- ggplot(g, aes(x=geneFamily, y=vntDensity, fill=acc)) + 
    xlab("gene familesi") + ylab("snp density per 10kb") +
    geom_boxplot(alpha=0.2) + ylim(c(0,200)) + 
    facet_wrap(~acc) +
    opts(axis.text.x = theme_text(angle = 85, hjust = 1, size = 7, colour = "grey50"));
  ggsave(p, filename = file.path(DIR_Stat,paste("snpDensity",".png",sep="")), width=8, height=5);
}
plotG_LR <- function(g) {
  for (i in 1:length(geneNames)) {
    gg=g[g$geneFamily == geneNames[i],];
    title = geneNames[i];
    coefs <- ddply(gg, .(acc), function(df) { 
      m <- lm(vntDensity ~ avgCov, data=df) 
      data.frame(a = coef(m)[1], b = coef(m)[2]) 
    });
  p <- ggplot(gg, aes(x=avgCov, y=vntDensity)) + 
    xlim(c(0,20)) + xlab("Average Coverage") + ylim(c(0,250)) + 
    #geom_point(alpha=0.5) +
    geom_abline(data=coefs, aes(intercept=a, slope=b, color=acc)) +
    opts(title=title);
  ggsave(p, filename = file.path(DIR_Stat,paste("snpDensity_",title,".png",sep="")),
    width=8, height=5);
  }
}
agg3 <- function(d) {
  #aggregate by geneFamily & acc & type
  q=aggregate(d$totCov,by=list(d$geneFamily,d$type,d$acc),FUN=sum);
  w=aggregate(d$totVnt,by=list(d$geneFamily,d$type,d$acc),FUN=sum);
  e=aggregate(d$calledVnt,by=list(d$geneFamily,d$type,d$acc),FUN=sum);
  r=aggregate(d$coveredVnt,by=list(d$geneFamily,d$type,d$acc),FUN=sum);  
  j=aggregate(d$length,by=list(d$geneFamily,d$type,d$acc),FUN=sum);  
  a=cbind(q,w[4],e[4],r[4],j[4]);
  al=length(a);
  a[al+1]=a[4]/a[5];
  a[al+1][is.na(a[al+1])] <- 0;
  a[al+2]=a[6]/((a[8]/a[5])*a[7]) * 10000;
  a[al+2][is.na(a[al+2])] <- 0;
  colnames(a)=c('geneFamily', 'type', "acc", 'totCov', 'totVnt', 
    'calledVnt', 'coveredVnt', 'length', 'avgCov', 'vntDensity');
}
plotA <- function(a, subGeneNames=0, suBACcs=0) {
  pa = 85; pw = 8; ph = 5; ps1 = 8; ps2 = 10;
  if (subGeneNames != 0 & suBACcs != 0) {
    a <- subset(a, geneFamily %in% subGeneNames & acc %in% suBACcs, drop=TRUE);
    a$geneFamily <- factor(a$geneFamily);
    a$acc <- factor(a$acc);
    pa = 30; pw = 6; ph = 3;
  }
  p <- ggplot(a, aes(acc, avgCov)) + 
    xlab("Accessions") + ylab("Average Coverage") +
    geom_bar(aes(fill=geneFamily), stat="identity", position="dodge") + 
    facet_grid(~ type)  + scale_fill_brewer(palette="Set1") +
    opts(title = "SNP Coverage Distribution") + 
    opts(axis.text.x = theme_text(angle = pa, hjust = 1, size = ps1, colour = "grey50"));
  ggsave(p, filename = file.path(dir,"avgCov.png"), width=pw, height=ph);
  p <- ggplot(a, aes(acc, vntDensity, fill=geneFamily)) + 
    xlab("Accessions") + ylab("SNP Density (10kb)") + labs(fill='Gene Family') +
    geom_bar(stat="identity", position="dodge") + 
    facet_grid(~ type)  + scale_fill_brewer(palette="Set1") +
    opts(axis.text.x = theme_text(angle = pa, hjust = 1, size = ps2, colour = "grey50"));
  ggsave(p, filename = file.path(DIR_Stat,"snpDensity.png"), width=pw, height=ph);
}

plotChr <- function(fName, suBACcs) {
  #plot snp coverage and density along chromosome
  xlabelSize = rep(10,8);
  d <- read.table(fName, header=TRUE);
  ld = length(d);
  d[ld+1]=d$totCov/d$totVnt;
  d[ld+1][is.na(d[ld+1])] <- 0;
  d[ld+2]=d$calledVnt/((d$length/d$totVnt) * d$coveredVnt) * 10000;
  #d[ld+2][is.na(d[ld+2])] <- 0;
  d$posMiddle = d$posMiddle / 100000;
  colnames(d)=c('chr', 'posMiddle', 'region', "acc", 'totCov', 
    'totVnt', 'calledVnt', 'coveredVnt', 'length', 'avgCov', 'vntDensity');
  d0 <- d;
  d <- subset(d, acc %in% suBACcs, drop=TRUE);
  d$acc <- factor(d$acc);
  cntAcc = 0;
  for (accs in levels(d$acc)) {
    cntAcc = cntAcc + 1;
    subD <- subset(d, d$acc==accs);
    subD <- subD[order(subD$posMiddle),];
    x <- c();
    for (j in 1:length(subD$vntDensity)) {
      idxStart = 1; idxEnd = length(subD$vntDensity);
      if (j-idxStart>15) {idxStart = j-15;}
      if (idxEnd-j>14) {idxEnd = j+14;}
      cav = sum(subD$calledVnt[idxStart:idxEnd]);
      le = sum(subD$length[idxStart:idxEnd]);
      cov = sum(subD$coveredVnt[idxStart:idxEnd]);
      tov = sum(subD$totVnt[idxStart:idxEnd]);
      vd = cav / ((le / tov) * cov) * 10000;
      x <- c(x, vd);
    }
    newd = data.frame(acc=accs,posMiddle=subD$posMiddle, vntDensity=x);
    if (cntAcc == 1) { e <- newd; }
    else { e <- rbind(e,newd); }
  }
  p <- ggplot(e, aes(posMiddle, vntDensity, colour=acc)) + 
    xlab("Location on chromosome (/100kb)") + ylab("SNP density (per 10kb)") +
    labs(colour='Accession') +
    geom_line(size=0.3) + scale_color_brewer(palette="Set1") +
    opts(title=chrName) +
    opts(axis.text.x = theme_text(size = xlabelSize[i], colour = "grey50"));
  ggsave(p, filename = file.path(DIR_Stat,paste("snpDensity_",chrName,".png",sep="")),
    width=10, height=4);
  p <- ggplot(e, aes(posMiddle, avgCov, colour=acc)) + 
    xlab("Location on chromosome (/100kb)") + ylab("Coverage") +
    geom_line(size=0.3) + ylim(c(0,80)) + scale_fill_brewer(palette="Set1") +
    facet_grid(acc ~ .)  + scale_color_brewer(palette="Set1") +
    opts(title=chrName) +
    opts(axis.text.x = theme_text(size = xlabelSize[i], colour = "grey50"));
  ggsave(p, filename = file.path(DIR_Stat,paste("avgCov_",chrName,".png",sep="")),
    width=10, height=4);
}
snpEffect <- function() {
  f01 = file.path(DIR_Misc2, 'genes', '05_vntEffect.txt');
  d1 <- read.table(f01, sep='\t', quote='"', header=TRUE);
  f02 = file.path(DIR_In, 'gene_id_r.txt');
  d2 <- read.table(f02, sep='\t', header=TRUE);
  f03 = file.path(DIR_In, 'go_info.txt');
  d3 <- read.table(f03, sep='\t', header=FALSE);
  colnames(d3) = c('family', 'desc');
  des = merge(d2, d1, by.x='gene', by.y='gene');
  des.1 <- cbind(des, 1); 
  colnames(des.1)[length(des.1)] <- "cnt";

  subEffects = c('Frameshift');
  suBACcs = c('HM017', 'HM018', 'HM022', 'HM029', 'HM030');
  subFamilies = c(); #c("GO:0043167", "GO:0000166", "GO:0016787", "GO:0001882", "NB-ARC", "NCR");
  groupEffect = 1;
  groupAcc = 1;

  des.2 = des.1;
  if(length(suBACcs) >= 1) {
    des.2 = subset(des.1, acc %in% suBACcs, drop=TRUE);
    des.2$acc <- factor(des.2$acc);
  }
  des.3 = des.2;
  if(length(subEffects) >= 1) {
    des.3 <- subset(des.2, ! effect %in% subEffects, drop=TRUE);
    des.3$effect <- factor(des.3$effect);
  }
  des.4 = des.3;
  if(length(subFamilies) >= 1) {
    des.4 <- subset(des.3, family %in% subFamilies, drop=TRUE);
    des.4$family <- factor(des.4$family);
  }
  
  des.5 = des.4;
  if(groupAcc == 1) {
    des.5$acc = NULL;
    des.5 = unique(des.5);
    des.6 = aggregate(des.5$cnt, by = list(des.5$family, des.5$effect), FUN=sum);
    colnames(des.6) = c('family', 'effect', 'cnt');
  } else {
    des.6 = aggregate(des.5$cnt, by = list(des.5$family, des.5$effect, des.5$acc), FUN=sum);
    colnames(des.6) = c('family', 'effect', 'acc', 'cnt');
  }
  
  des.7 = des.6;
  if(groupEffect == 1) {
    effectH = data.frame(cbind(effectD = c("3' UTR", "5' UTR", "Essential splice site", "Intronic", "Non-synonymous", "Splice site", "Start lost", "Stop gained", "Stop lost", "Synonymous"), effectS =  c("UTR", "UTR", "Large-effect", "Intronic", "Non-synonymous", "Non-synonymous", "Large-effect", "Large-effect", "Large-effect", "Synonymous")));
    des.71 <- merge(des.6, effectH, by.x='effect', by.y='effectD');
    des.7 <- aggregate(des.71$cnt, by = list(des.71$family, des.71$effectS), FUN=sum);
    colnames(des.7) = c('family', 'effect', 'cnt');
  }

  des.8 <- cbind(des.7, 0);
  tmpLen <- length(colnames(des.8));
  colnames(des.8)[tmpLen] <- c("percent");
  for (cc in 1:length(des.8[,1])) {
    if(groupAcc == 1) {
      sumT = sum(des.8$cnt[des.8$family == as.character(des.8$family[cc])]);
    } else {
      sumT = sum(des.8$cnt[des.8$family == as.character(des.8$family[cc]) & des.8$acc==as.character(des.8$acc[cc])]);
    }
    des.8$percent[cc] = des.8$cnt[cc] / sumT;
  }
  if(length(suBACcs) == 0 & groupAcc == 0) {
    numcol = 4; yfontsize = 6; iw = 20; ih = 20;
  } else{
    numcol = length(suBACcs); yfontsize = 6; iw = 8; ih = 4;
  }
  lorder = c('Non-synonymous', 'Large-effect', 'Synonymous', 'Intronic', 'UTR');
  p <- ggplot(des.8, aes(family, percent, fill=factor(effect, levels=lorder, ordered=T))) + 
    xlab("Gene Family") + ylab("Proportion") + labs(fill='SNP effects') +
    geom_bar() + 
    coord_flip() +
    #facet_wrap( ~ acc, ncol=numcol ) +
    scale_fill_brewer(palette='Set1') +
    opts(title="SNP effect distribution [outgroup]") +
    opts(axis.text.y=theme_text(size=yfontsize, hjust=1, colour="royalblue"));
  ggsave(p, filename = file.path(DIR_Stat,paste("snpEffects_outgroup",".png",sep="")),
    width=iw, height=ih);
}
plotGeneStat <- function() {
  f1 = file.path(DIR_Misc1, "cat", "31_cat_info.txt")
  f2 = file.path(DIR_Misc1, "cat", "32_cat.txt")
  f3 = file.path(DIR_Misc1, "cat", "51_stat.txt")
  catD = read.table(f1, sep="\t", header=TRUE)
  cat = read.table(f2, sep="\t", header=TRUE)
  stat = read.table(f3, sep="\t", header=TRUE)
  stat = merge(stat, cat, by.x='id', by.y='id')
  stat = merge(stat, catD[,c(1,3)], by.x='family', by.y='familyId')

  famIdsDel = catD$familyId[grep("(reverse transcriptase)|(IMP dehydrogenase)|(Ribonuclease)|(Polynucleotidyl transferase)", catD$family)]
  stat = stat[! stat$family %in% famIdsDel, ]

  col = "ThetaPi"
  statM = aggregate(stat[,col], by=list(stat$family), median, na.rm=TRUE)
  colnames(statM) = c("familyId", col)
  catM = merge(catD, statM, by.x="familyId", by.y="familyId")
  catStr = paste(catM$family, " (n=", catM$nGene, ")", sep="");
  catM = cbind(catM, catStr=catStr)
  catM = catM[order(-catM[,col]), ]
  stat$family = factor(stat$family, levels=catM$familyId)
  stat$dataset = factor(stat$dataset, levels=sort(unique(stat$dataset)))
  p <- ggplot(stat) + 
    geom_boxplot(aes(family, stat[,col], group=family, fill=dataset)) +
    scale_x_discrete(name='', breaks=catM$familyId, labels=catM$catStr) +
    scale_y_continuous(name='', limits=c(0,0.03)) +
    scale_fill_manual(breaks=seq(1,4,1), labels = c('Gene Ontology','Domain Function', 'NBS-LRRs', 'CRPs'), values=c('yellow', 'orange', 'royalblue', 'red')) +
    coord_flip() +
    opts(legend.title = theme_blank()) +
    opts(axis.text.y = theme_text(size=8, hjust=1, colour='blue'));
  ggsave(p, filename=file.path(DIR_R, 'cat', paste(col, ".png", sep="")), width=10, height=7);
}
plotRsq <- function(fName) {
  #plot r2 with distance
  r <- read.table(file.path(dir, fName), header=TRUE);
  #rsq = as.numeric(levels(r$rsq)[as.integer(r$rsq)]);
  #distance = as.numeric(levels(r$distance)[as.integer(r$distance)]);
  #r <- as.data.frame(cbind(rsq,distance));
  fit <- lm(rsq~distance,data=r);
  plot(r$distance,r$rsq);
  abline(fit);
}
plotVnt <- function() {
  f01 = file.path(DIR_Misc2, 'genes', '06_vntCnt.txt');
  d <- read.table(f01, sep = "\t", header = TRUE);
  d <- cbind(d, 0, 1);
  colnames(d)[1] = "family";
  colnames(d)[(length(d)-1):length(d)] = c('vd', 'cnt');
  d$vd = d$vntCnt / (d$len * d$pctUniqCov);

  d.s = d;
  suBACcs = c(); #c('HM005', 'HM006', 'HM015', 'HM101');
  subFamilies = c(); #c("GO:0043167", "GO:0000166", "GO:0016787", "GO:0001882", "NB-ARC", "NCR");
  if(length(suBACcs) >= 1) {
    d.s = subset(d.s, acc %in% suBACcs, drop=TRUE);
    d.s$acc <- factor(d.s$acc);
  }
  if(length(subFamilies) >= 1) {
    d.s <- subset(d.s, family %in% subFamilies, drop=TRUE);
    d.s$family <- factor(d.s$family);
  }
  tmp1 = aggregate(d.s$avgCov, by=list(d.s$family, d.s$acc), FUN=mean, na.rm=TRUE);
  tmp2 = aggregate(d.s$avgCov, by=list(d.s$family, d.s$acc), FUN=sd, na.rm=TRUE);
  tmp3 = aggregate(d.s$vd, by=list(d.s$family, d.s$acc), FUN=mean, na.rm=TRUE);
  tmp4 = aggregate(d.s$vd, by=list(d.s$family, d.s$acc), FUN=sd, na.rm=TRUE);
  tmp5 = aggregate(d.s$cnt, by=list(d.s$family, d.s$acc), FUN=sum);
  stat = cbind(tmp1, tmp2[ ,3], tmp3[ ,3], tmp4[ ,3], tmp5[ ,3], 0, 0);
  colnames(stat) =   c('family', 'acc', 'covMean', 'covSd', 'vdMean', 'vdSd', 
    'cnt', 'covSe', 'vdSe');
  stat$covSe = stat$covSd * 1.96 / sqrt(stat$cnt);
  stat$vdSe = stat$vdSd * 1.96 / sqrt(stat$cnt);
  limitsCov <- aes(ymax = covMean + covSe, ymin = covMean - covSe);
  limitsVd <- aes(ymax = vdMean + vdSe, ymin = vdMean - vdSe);
  p <- ggplot(stat, aes(acc, covMean, color = family)) + 
    xlab("accession") + ylab("uniq coverage") + 
    geom_bar(position='dodge', fill = 'white') + 
    geom_errorbar(limitsCov, position='dodge', width=0.9);
  ggsave(p, file = file.path(DIR_Stat, "cov.png"), width=13, height = 7);
  p <- ggplot(stat, aes(acc, vdMean, color = family)) + 
    xlab("accession") + ylab("variant density") +
    geom_bar(position='dodge', fill = 'white') + 
    geom_errorbar(limitsVd, position='dodge', width=0.9);
  ggsave(p, file = file.path(DIR_Stat, "vd.png"), width=13, height = 7);
  if(length(suBACcs) == 0) {
    numcol = 4; ymax = 50; xfontsize = 6; xfontangle = 30; xhjust = 1; iw = 20; ih = 20;
  } else{
    numcol = length(suBACcs); ymax = 40; xfontsize = 6; xfontangle = 30; xhjust = 1; 
    iw = 6; ih = 4;
  }
  p <- ggplot(d.s, aes(family, avgCov, fill = family)) + 
    ylab("coverage (fold)") + ylim(c(0, ymax)) + xlab('') +
    geom_boxplot() + facet_wrap( ~ acc, ncol=numcol) +
    opts(title = "Sequencing Depth") +
    opts(axis.text.x=theme_text(angle=xfontangle, hjust=xhjust, size=xfontsize, colour="blue"));
  ggsave(p, file = file.path(DIR_Stat, "avgCov.png"), width=iw, height = ih);
  p <- ggplot(d.s, aes(family, avgCov, fill = family)) + 
    ylab("Uniq Coverage (fold)") + ylim(c(0, ymax)) + xlab('') +
    opts(title = "Sequencing Depth") +
    geom_boxplot() + facet_wrap( ~ acc, ncol=numcol) +
    opts(axis.text.x=theme_text(angle=xfontangle, hjust=xhjust, size=xfontsize, colour="blue"));
  ggsave(p, file = file.path(DIR_Stat, "avgUniqCov.png"), width=iw, height = ih);
  p <- ggplot(d.s, aes(family, vd, fill = family)) + 
    ylab("Variant Density (/bp)") + ylim(c(0, 0.04)) + xlab('') +
    opts(title = "Variant Density distribution") +
    geom_boxplot() + facet_wrap( ~ acc, ncol=numcol) +
    opts(axis.text.x=theme_text(angle=xfontangle, hjust=xhjust, size=xfontsize, colour="blue"));
  ggsave(p, file = file.path(DIR_Stat, "vd_uniqCov.png"), width=iw, height = ih);
}
plotDensity1 <- function() {
  dirI = file.path(DIR_Misc1, "circos", "06_density")
  den1 = read.table(file.path(dirI, "01_ncr.txt"), sep="\t", header=FALSE)
  colnames(den1) = c("chr", "start", "end", "count")
  den2 = read.table(file.path(dirI, "02_te.txt"), sep="\t", header=FALSE)
  colnames(den2) = c("chr", "start", "end", "count")
  p <- ggplot(data=den1) + 
    geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=count, fill="CRP"), size=0, alpha=1) +
    layer(data=den2, geom='rect', mapping=aes(xmin=start, xmax=end, ymin=0, ymax=count, fill="TE"), geom_params=list(size=0, alpha=0.3)) +
    layer(data=a, geom='rect', mapping=aes(xmin=start, xmax=end, ymin=-1, ymax=0, fill=type), geom_params=list(size=0)) +
    scale_fill_manual(values=c('CRP'='forestgreen', 'TE'='royalblue', 'phase 1 BAC'='orchid1', 'phase 2 BAC'='salmon', 'phase 3 BAC'='red1', 'centromere'='black', 'BAC of unknown phase'='green'), legend=TRUE) +
    scale_colour_manual(values=c('firebrick', 'forestgreen')) +
    scale_x_continuous(name='chr position (bp)', formatter='comma') +
    scale_y_continuous() +
    facet_grid(chr ~ .) + 
    opts(legend.title=theme_blank());
  ggsave(p, filename = file.path(DIR_R, 'crp.png'), width=12, height=8);
}



summarise.interval=function(rates.file="rates.txt", burn.in=30, locs.file=FALSE) {
  x = read.table(rates.file, skip=1, fill=T);
  x = as.matrix(x);
  low = as.integer(nrow(x)*burn.in/100);
#  cat("\n\nSummarise output from MCMC estimation of recombination rates in INTERVAL (LDhat 2.1)\n\n");
  cat(paste("Number of SNPs = ", ncol(x), "\n", sep=""));
  cat(paste("Number of samples = ", nrow(x), "\n", sep=""));
  cat(paste("Burn-in period = ", low, " samples\n", sep=""));
#  x11();
#  par(mfrow=c(1,2));
#  plot(x[,1], type="s", col=rgb(0,0,0.5), xlab="Sample", ylab="Total map length", 
#  main="Mixing of total map length");
#  image(x=c(1:nrow(x)), y=c(1:(ncol(x)-1)), z=log(x[,2:ncol(x)]), xlab="Sample", 
#  ylab="log(rate) at SNP", main="Mixing of rates");

  means<-apply(x[low:nrow(x),], 2, mean, na.rm=T);
  q.95<-apply(x[low:nrow(x),], 2, quantile, probs=c(0.025, 0.5, 0.975), na.rm=T);
  if (locs.file==FALSE) {pos<-c(1:ncol(x)); xlab<-"Position (SNP)";}
  else {pos<-as.vector(as.matrix(read.table(locs.file, as.is=T, skip=1))); xlab<-"Position (Relative, /kb)"}
  tmp <- 0;
  for(j in 2:length(means)) {
    tmp <- tmp + means[j]*(pos[j]-pos[j-1]);
  }
  cat(paste("Mean posterior total map length (4Ner) = ", signif(means[1], 4), "(", signif(tmp,4), ")\n\n", sep=""));
  
#  x11();
#  plot(pos[1:(length(pos)-1)], y=means[2:length(means)], type="s", col=rgb(0,0,0.5),
#  xlab=xlab, ylab="Posterior mean rate", main="Posterior mean rates");
#  lines(pos[1:(length(pos)-1)], y=q.95[1,2:length(means)], type="s", col=grey(0.75), lty="dotted");
#  lines(pos[1:(length(pos)-1)], y=q.95[3,2:length(means)], type="s", col=grey(0.75), lty="dotted");
  op<-cbind(means, t(q.95));
  colnames(op)<-c("Mean", "q2.5", "Median", "q97.5");
  return(op);
}
summarise.rhomap<-function(rates.file="rates.txt", burn.in=30, locs.file=FALSE) {
  return(summarise.interval(rates.file=rates.file, burn.in=burn.in, locs.file=locs.file));
}
summarise.pairwise<-function(surf=TRUE, window=FALSE, rm=FALSE, test=FALSE, ci=FALSE, locs.file=FALSE) {
  cat("\n\nSummarising output from PAIRWISE program in LDhat 2.1\n\n");
  if (locs.file==TRUE) pos<-as.vector(as.matrix(read.table(locs.file, as.is=T, skip=1)));
  if (surf) {
    surface = read.table("outfile.txt", skip=10, fill=T);
    x11();
    plot(surface[,1:2], type="l", col=rgb(0,0,0.5), xlab="4Ner", 
    ylab="Composite likelihood", main="Likelihood surface for 4Ner over region");
    xmax=surface[1,1];
    ymax = surface[1,2];
    for (i in 2:nrow(surface)) {
      if (surface[i,2]>ymax) {ymax=surface[i,2]; xmax=surface[i,1];}
    }
    cat(paste("Maximum composite likelihood estimate of 4Ner = ", xmax, "\n", sep=""));
  }
  if (window) {
    win = read.table("window_out.txt", skip=5);
    x11();
    par(mfrow=c(3,1));
    plot(x=(win[,1]+win[,2])/2, win[,3], type="n", xlab="Window", ylab="Number SNPs", main="SNP density");
    segments(x0=win[,1], y0=win[,3], x1=win[,2], y1=win[,3], col="blue");
    plot(x=(win[,1]+win[,2])/2, win[,4], type="n", xlab="Window", ylab="4Ner per kb/bp", main="Local recombination rate");
    segments(x0=win[,1], y0=win[,4], x1=win[,2], y1=win[,4], col="red");

    av = vector(length=nrow(win));
    for(i in 1:nrow(win)) av[i]=xmax/(win[nrow(win),2]-win[1,1]);
    lines(x=(win[,1]+win[,2])/2, y=av, col="green");
    plot(x=(win[,1]+win[,2])/2, y=win[,5], type="n", xlab="Window midpoint", ylab="CLR", main="Composite likelihood ratio");
    segments(x0=win[,1], y0=win[,5], x1=win[,2], y1=win[,5], col="black");
  }
  if (rm) {
    rmin = read.table("rmin.txt", skip=5, fill=T);
    rmin= as.matrix(rmin[1:nrow(rmin),2:ncol(rmin)]);

    f<-as.matrix(read.table("freqs.txt", as.is=T, skip=5)[,2:6]);
    n.chr<-sum(f[1,]);
    anc<-apply(f[,2:5], 1, sum);
    mx<-apply(f[,2:5], 1, max);
    seg<-mx<anc;
    
    rmin<-rmin[seg[1:(length(seg)-1)],];
    rmin<-rmin[,seg[2:length(seg)]];

    maxrm=0;
    maxscrm=0;
    for (i in 1:nrow(rmin)) {
      if (i<ncol(rmin)) {for (j in (i):ncol(rmin)) {if (rmin[i,j]>maxrm) maxrm=rmin[i,j];}}
      if (i>1) {for (j in 1:(i-1)) {if(rmin[i,j]>maxscrm) maxscrm=rmin[i,j];}}
    }

    xu=matrix(0,nrow=nrow(rmin)+1, ncol=nrow(rmin)+1);
    xl=matrix(0,nrow=nrow(rmin)+1, ncol=nrow(rmin)+1);
    
    for (i in 1:nrow(rmin)) {
      if (i<ncol(rmin)) {
        for (j in (i):ncol(rmin)) {xu[i,j+1]=rmin[i,j]; xl[i,j+1]=0;}
      }
      if (i>1) {
        for (j in 1:(i-1)) {xl[i,j]=rmin[i,j]; xu[i,j+1]=0;}
      }
      xu[i,i]=0;
      xl[i,i]=0;
    }

    blues = c(rgb(0,0,(15:6)/15));
    reds = c(rgb((6:15)/15,0,0));

    x11();
    par(mfrow=c(2,1));
    
    image(xu,  xaxt="n", yaxt="n", xlab="Position", ylab="Position", main="Total recombination", 
  col=c("white", blues, reds, "yellow"));
    image(log(xl+0.001), main="Recombination scaled by distance", xaxt="n", yaxt="n", 
  xlab="Position", ylab="Position", col=c("white", blues, reds, "yellow"), 
  breaks=as.vector(quantile(xl, probs=seq(0,1,length.out=23))));
    rm(xu);
    rm(xl);
    rm(rmin);
  }

  if (test) {
    xx<-max(surface[,2]);
    rd<-as.matrix(read.table("rdist.txt", skip=4)[,2:5]);

    x11();
    hist(rd[,2], col="blue", breaks=seq(xlim[1], xlim[2], length.out=20), , 
  main="Likelihood permutation test", xlab="Composite Likelihood");
    points(xx,1,pch=25, col="red", bg="red");
  }

  if (ci) {
  cat("\n\nNothing implemented yet to display sampling distribution\n\n");
  }
}
