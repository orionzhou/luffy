
### function to cluster a vector of sorted positions using a given window size (no need to modify)
locCluster <- function(pos, wsize) {
  npos = length(pos)
  df = data.frame(id=1:npos, pos=pos, cluster=1:npos)
  df = df[order(df$pos),]
  for (i in 1:npos) {
    for (j in (i+1):npos) {
      if(j>npos) {
        next
      }
      if(df$pos[j] - df$pos[i] <= wsize) {
        df$cluster[j] = df$cluster[i]
      } else {
        next
      }
    }
  }
  clusters = unique(df$cluster)
  tmp = cbind(cluster1=clusters, cluster2=1:length(clusters))
  x = merge(df, tmp, by.x='cluster', by.y='cluster1')
  
  df2 = data.frame(id=x$id, pos=x$pos, cluster=x$cluster2, cluster_y=1)
  df2 = df2[order(df2$pos),]
  
  clusterP = ''
  for (i in 1:npos) {
    if(df2$cluster[i] == clusterP) {
      df2$cluster_y[i] = df2$cluster_y[i-1] + 1
    } else {
      clusterP = df2$cluster[i]
    }
  }

#  hist(as.matrix(table(df2$cluster)), xlab='cluster size', main=paste(wsize, 'bp', sep=''))
  df2
}

### create a data frame of genomic locations
tt = data.frame(
	chr = c('chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr2'),
	beg = c(100, 115, 150, 170, 190, 180, 250, 271),
	end = c(100, 115, 150, 170, 190, 180, 250, 271)
)

dirw = '/home/springer/zhoux379/data/in'
fi = file.path(dirw, 'test1.tsv')
ti = read.table(ft, header = T, sep = "\t", as.is = T)
tt = ti[,3:5]

### compute a global position for each row - assume each chromosome is 100Mb long so no inter-chromosomal locations get clustered
chrs = sort(unique(tt$chr))
chrmap = 1:length(chrs)
names(chrmap) = chrs

tt = cbind(tt, pos = chrmap[tt$chr] * 100000000 + (tt$beg + tt$end) / 2)
tt = tt[order(tt$pos), ]

### This is what "tt" looks like
#    chr beg end       pos
# 1 chr1 100 100 100000100
# 2 chr1 115 115 100000115
# 3 chr1 150 150 100000150
# 4 chr1 170 170 100000170
# 5 chr1 190 190 100000190
# 6 chr2 180 180 200000180
# 7 chr2 250 250 200000250
# 8 chr2 271 271 200000271

tt2 = locCluster(tt$pos, 200000)

### output - "cluster": cluster ID; "cluster_y": numeric ID within one cluster
#   id       pos cluster cluster_y
# 1  1 100000100       1         1
# 2  2 100000115       1         2
# 3  3 100000150       2         1
# 4  4 100000170       2         2
# 5  5 100000190       2         3
# 6  6 200000180       3         1
# 7  7 200000250       4         1
# 8  8 200000271       5         1

fo = file.path(dirw, 'test1.out.tsv')
write.table(tt2[,3:4], fo, sep = "\t", row.names = F, col.names = T, quote = F)
