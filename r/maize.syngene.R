require(tidyverse)

dirw = '/home/springer/zhoux379/data/genome/B73/gene_mapping'
fi = file.path(dirw, "maize_syntenic_genes_Alex\'s_paper.txt")
ti = read_tsv(fi)

ti = ti[,1:12]
colnames(ti) = c("achrom", "abeg", "aend", "agid", "bchrom1", "bbeg1", "bend1", "bgid1", "bchrom2", "bbeg2", "bend2", "bgid2")

gids1 = ti$bgid1[!is.na(ti$bgid1)]
gids2 = ti$bgid2[!is.na(ti$bgid2)]
gids = unique(c(gids1, gids2))
stopifnot(length(gids) == length(unique(gids1)) + length(unique(gids2)))

to1 = tibble(gid = gids1, subgenome = 'subgenome1')
to2 = tibble(gid = gids2, subgenome = 'subgenome2')
to = rbind(to1, to2)

fo = file.path(dirw, "syn.gid.tsv")
write_tsv(to, fo)
