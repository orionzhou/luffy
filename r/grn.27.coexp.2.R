source("br.fun.R")

#dirw = file.path(Sys.getenv("misc2"), "grn23", "47.coexp.test")
dirw = file.path(Sys.getenv("misc2"), "briggs", "47.coexp.test")

#fi = file.path(dirw, "../37.rpkm.filtered.tsv")
fi = file.path(dirw, "../36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)
ti = spread(ti[ti$Genotype == "B73", c(1,2,5)], Tissue, fpkm)

expr = t(as.matrix(ti[,-1]))
colnames(expr) = ti[,1]
gids = ti[,1]
ng = length(gids)



fc = file.path(dirw, "params.tsv")
tc = read.table(fc, sep = "\t", header = T, as.is = T)
tce = unique(tc[,2:3])

dummy <- matrix(rep(NA, ng*ng), nrow=ng)
ii = which(lower.tri(dummy))
colidx = as.integer((ii - 1) / ng) + 1
rowidx = ii - (colidx - 1) * ng

for (i in 1:nrow(tce)) {
	coex_opt = tce$coex_opt[i]; n_edges = tce$n_edges[i]	
	net = sprintf("%s_%dM", coex_opt, n_edges)

	fd = sprintf("%s/01.edgeweight/%s.rda", dirw, coex_opt)
	if(i %in% c(1, 3, 5)) {
		x = load(fd)
		if(coex_opt == 'R') {
			ord = order(coexv)
		} else {
			ord = order(-coexv)
		}
	}
	idxs = ord[1:(n_edges*1000000)]
	dw = data.frame(g1 = gids[rowidx[idxs]], g2 = gids[colidx[idxs]], coex = coexv[idxs])
	#length(unique(c(dw$g1, dw$g2)))
	if(coex_opt == 'R') {
		dw$coex = exp((-(dw$coex-1)/25))
	}
	
	fo = sprintf("%s/02.edges/%s.tsv", dirw, net)
	write.table(dw, fo, sep = "\t", row.names = F, col.names = F, quote = F)
}

to = data.frame()
for (i in 1:nrow(tce)) {
	coex_opt = tce$coex_opt[i]; n_edges = tce$n_edges[i]	
	net = sprintf("%s_%dM", coex_opt, n_edges)

	fd = sprintf("%s/02.edges/%s.tsv", dirw, net)
	td = read.table(fd, header = F, as.is = T, sep = "\t")
	cat(sprintf("%s: %d unique genes\n", net, length(unique(c(td$V1, td$V2)))))
	
	g = as.undirected(graph(c(t(as.matrix(td[,1:2])))))
	#E(g)$weight = td$V3
	V(g)$size <- 2
	V(g)$frame.color <- "white"
	V(g)$color <- "orange"
	V(g)$label <- "" 
	E(g)$arrow.mode <- 0
	
	fp = sprintf("%s/%s.pdf", dirw, coex_opt)
	#plot(g, layout = layout_in_circle)
	deg = degree(g)
	td1 = as.data.frame(table(deg))
	colnames(td1) = c("degree", "freq")
	td1$degree = as.numeric(levels(td1$degree))[td1$degree]
	td1 = cbind(net = net, td1)
	to = rbind(to, td1)
}
tos = ddply(to, .(net), summarise, ngene = sum(freq))
netnames = sprintf("%s: %d genes", tos$net, tos$ngene)
names(netnames) = tos$net

fp = sprintf("%s/02.edges/degree_dist.pdf", dirw)
p1 = ggplot(to, aes(x = degree, y = freq)) +
  	geom_line(stat = 'identity') +
	#scale_x_continuous(name = '') +
	#scale_y_continuous(expand = c(0.08, 0)) +
	#scale_fill_manual(values = c("firebrick", "grey"), labels = c("in modules", "not in modules"), name = net) +
	facet_wrap(~net, scale = 'free', nrow = 3, labeller = as_labeller(netnames)) +
	theme_bw() +
	#theme(legend.position = c(0.55, 0.04), legend.direction = "vertical", legend.justification = c(0,0), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 9), legend.box.background = element_rect(color='black')) +
	theme(axis.ticks.length = unit(0, 'lines')) +
	theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
	theme(axis.title.x = element_text(size = 9)) +
	theme(axis.title.y = element_text(size = 9)) +
	theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
	theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
ggsave(p1, filename = fp, width = 8, height = 8)