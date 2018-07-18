dirW = file.path(DIR_In, "csci8980/03_mdr_data")
dirO = file.path(dirW, "50_stat")

get_sum_stat = function(x) {
  x2 = unique(x[ , c('maf', 'h', 'size')])
  x2 = cbind(x2, time=NA, time_hi=NA, time_lo=NA, power1=NA, power1_hi=NA, power1_lo=NA, power2=NA, power2_hi=NA, power2_lo=NA)
  for (i in 1:nrow(x2)) {
    x.sub = subset(x, maf==x2$maf[i] & h==x2$h[i] & size==x2$size[i])
    
    time.mean = mean(x.sub$time)
    time.sd = sd(x.sub$time)
    x2$time[i] = time.mean
    x2$time_lo[i] = max(0, time.mean-time.sd)
    x2$time_hi[i] = time.mean+time.sd
    
    power1.mean = mean(x.sub$power1)
    power1.sd = sd(x.sub$power1)
    x2$power1[i] = power1.mean
    x2$power1_lo[i] = max(0, power1.mean-power1.sd)
    x2$power1_hi[i] = min(100, power1.mean+power1.sd)
    
    power2.mean = mean(x.sub$power2)
    power2.sd = sd(x.sub$power2)
    x2$power2[i] = power2.mean
    x2$power2_lo[i] = max(0, power2.mean-power2.sd)
    x2$power2_hi[i] = min(100, power2.mean+power2.sd)
  }
  x2$size=factor(x2$size)
  x2$h = factor(x2$h)
  x2$maf = factor(x2$maf)
  x2
}

f1 = file.path(dirW, "21_mdr.tbl")
a01 = read.table(f1, header=TRUE, sep="\t")
a11 = get_sum_stat(a01)

f2 = file.path(dirW, "22_plink.tbl")
a02 = read.table(f2, header=TRUE, sep="\t")
a12 = get_sum_stat(a02)

f4 = file.path(dirW, "24_team.tbl")
a04 = read.table(f4, header=TRUE, sep="\t")
a14 = get_sum_stat(a04)

f5 = file.path(dirW, "25_boost.tbl")
a05 = read.table(f5, header=TRUE, sep="\t")
a15 = get_sum_stat(a05)

plot_stat <- function(d, name) {
  p = ggplot(d) + 
    geom_bar(mapping = aes(x=size, y=time, fill=maf), stat='identity', position='dodge') + 
    geom_errorbar(mapping = aes(x=size, ymin=time_lo, ymax=time_hi), size=0.3, width=0.3, stat='identity', position=position_dodge(width=1)) + 
    facet_wrap( ~ h, nrow = 1) +
    scale_fill_brewer(palette='Set1') +
    scale_x_discrete(name="sample size") +
    scale_y_continuous(name="time(s)") +
    opts(title=name, axis.text.x=theme_text(size=8, angle=30));
  ggsave(p, filename=file.path(dirO, paste(name, "time.png", sep="_")), width=8, height=2.5)
  
  p = ggplot(d) + 
    geom_bar(mapping = aes(x=size, y=power1, fill=maf), stat='identity', position='dodge') + 
    geom_errorbar(mapping = aes(x=size, ymin=power1_lo, ymax=power1_hi), size=0.3, width=0.3, stat='identity', position=position_dodge(width=1)) + 
    facet_wrap( ~ h, nrow = 1) +
    scale_fill_brewer(palette='Set1') +
    scale_x_discrete(name="sample size") +
    scale_y_continuous(name="Power") +
    opts(title=name, axis.text.x=theme_text(size=8, angle=30));
  ggsave(p, filename=file.path(dirO, paste(name, "power1.png", sep="_")), width=8, height=2.5)

  p = ggplot(d) + 
    geom_bar(mapping = aes(x=size, y=power2, fill=maf), stat='identity', position='dodge') + 
    geom_errorbar(mapping = aes(x=size, ymin=power2_lo, ymax=power2_hi), size=0.3, width=0.3, stat='identity', position=position_dodge(width=1)) + 
    facet_wrap( ~ h, nrow = 1) +
    scale_fill_brewer(palette='Set1') +
    scale_x_discrete(name="sample size") +
    scale_y_continuous(name="Extended Power") +
    opts(title=name, axis.text.x=theme_text(size=8, angle=30));
  ggsave(p, filename=file.path(dirO, paste(name, "power2.png", sep="_")), width=8, height=2.5)
}

plot_stat(a11, "MDR")
plot_stat(a12, "PLINK")
plot_stat(a14, "TEAM")
plot_stat(a15, "BOOST")
