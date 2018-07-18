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

for (i in 1:24) {

	nid = tc$nid[i]; coex_opt = tc$coex_opt[i]; n_edges = tc$n_edges[i]
	algorithm = tc$algorithm[i]; mcl_param = tc$mcl_param[i]
	
	net = sprintf("%s_%dM_%s", coex_opt, n_edges, algorithm)
	if(algorithm == "mcl") {
		net = sprintf("%s_%g", net, mcl_param)
	}
	
	f_edge = sprintf("%s/02.edges/%s_%dM.tsv", dirw, coex_opt, n_edges)
	if(algorithm == 'mcl') {
		fo2 = sprintf("%s/05.modules/%s.2.txt", dirw, net)
		cmd = sprintf("mcl %s --abc -scheme 7 -I %g -o %s", f_edge, mcl_param, fo2)
		system(cmd)
	
		fo3 = sprintf("%s/05.modules/%s.tsv", dirw, net)
		cmd = sprintf("mcl2tsv.py %s %s", fo2, fo3)
		system(cmd)
	} else {
		stopifnot(algorithm == 'clusterone')
		fo2 = sprintf("%s/05.modules/%s.2.csv", dirw, net)
		cmd = sprintf("java -jar $src/cluster_one-1.0.jar -f edge_list -F csv %s > %s", f_edge, fo2)
		system(cmd)

		fo3 = sprintf("%s/05.modules/%s.tsv", dirw, net)
		cmd = sprintf("one2tsv.py %s %s", fo2, fo3)
		system(cmd)
	}
	
	fm = sprintf("%s/05.modules/%s.tsv", dirw, net)
	tm = read.table(fm, header = F, as.is = T, sep = "\t")
	cat(sprintf("%s: %d modules with %d genes\n", net, length(unique(tm$V1)), nrow(tm)))
}

for (i in 1:nrow(tc)) {
	nid = tc$nid[i]; coex_opt = tc$coex_opt[i]; n_edges = tc$n_edges[i]
	algorithm = tc$algorithm[i]; mcl_param = tc$mcl_param[i]
	
	net = sprintf("%s_%dM_%s", coex_opt, n_edges, algorithm)
	if(algorithm == "mcl") {
		net = sprintf("%s_%g", net, mcl_param)
	}
	fm = sprintf("%s/05.modules/%s.tsv", dirw, net)
	tm = read.table(fm, header = F, as.is = T, sep = "\t")
	cat(sprintf("%s: %d modules with %d genes\n", net, length(unique(tm$V1)), nrow(tm)))
}