require(dplyr)
require(ggmap)
require(xlsx)

dirw = file.path(Sys.getenv("misc2"), "urt")

### extract phenotypic data
fi = file.path(dirw, "00.urt.rda")
x = load(fi)

tnames = names(ud)
yrs = sapply(strsplit(tnames, split="[_]"), "[", 1)
exs = sapply(strsplit(tnames, split="[_]"), "[", 2)
phs = sapply(strsplit(tnames, split="[_]"), "[", 3)

tp = data.frame()
phs1 = c("Height", "Lodging", "Maturity", "Oil", "Protein", "SeedQuality", "SeedSize", "YieldBuA", "YieldRank")
phs2 = c("2YearMean", "3YearMean", "4YearMean", "RegionalSummary", "Entries")
phs3 = c("Descriptive", "DescriptiveDisease")
for (i in 1:length(ud)) {
	yr = yrs[i]; ex = exs[i]; ph = phs[i]
	#DescriptiveDisease need separate processing
	if(ph %in% phs1) {
		tp1 = ud[[i]]
		#if("Mean" %in% colnames(tp1)) { tp1 = tp1[,-which(colnames(tp1)=="Mean")] }
		#if(ph == "Maturity") { for (k in 2:ncol(tp1)) { tp1[,k] = as.character(tp1[,k]) }}
		if(colnames(tp1)[1] != 'Strain') {
			cat(tnames[i], "not start with Strain: ", colnames(tp1)[1], "\n")
			colnames(tp1)[1] = "Strain"
		}
		tp2 = reshape(tp1, direction='long', varying=list(2:ncol(tp1)), idvar='Strain', timevar="Location", v.names='Value', times=colnames(tp1)[2:ncol(tp1)])
		tp3 = cbind.data.frame(Year=yr, Test=ex, Phenotype=ph, tp2, stringsAsFactors=F)
		tp3$Value = as.character(tp3$Value)
		tp = rbind(tp, tp3)
	} else if(ph %in% phs2) {
	} else if(ph %in% phs3) {
		tp1 = ud[[i]]
		if(colnames(tp1)[1] != 'Strain') {
			cat(tnames[i], "not start with Strain: ", colnames(tp1)[1], "\n")
			colnames(tp1)[1] = "Strain"
		}
		if(colnames(tp1)[2] != 'DescriptiveCode') {
			cat(tnames[i], "not start with DescriptiveCode: ", colnames(tp1)[2], "\n")
		} else {
			ph = 'DescriptiveCode'
			tp3 = data.frame(Year=yr, Test=ex, Phenotype=ph, Strain = tp1$Strain, Location='Mean', Value=as.character(tp1[,2]), stringsAsFactors = F)
			tp = rbind(tp, tp3)
		}
		if(ncol(tp1) <= 2) next
		zz = strsplit(colnames(tp1)[3:ncol(tp1)], split="[_]")
		for (j in 1:length(zz)) {
			ary = zz[[j]]
			if(length(ary) < 2) next
			if(length(ary) > 2) {
				cat(tnames[i], " col error: ", colnames(tp1)[j+2], "\n")
				if(ary[2] == '' | ary[2] == ary[3]) {
					ary[2] = ary[3]
				} else {
					stop("fatal error: unknown col")
				}
			}
			ph = ary[1]; lc = ary[2]
			tp3 = data.frame(Year=yr, Test=ex, Phenotype=ph, Strain = tp1$Strain, Location=lc, Value=as.character(tp1[,2+j]), stringsAsFactors = F)
			tp = rbind(tp, tp3)
		}
	} else {
		cat("unknown phenotype: ", ph, "\n")
	}
}

fo = file.path(dirw, "01.tsv")
write.table(tp, fo, sep = "\t", row.names = F, col.names = T, quote = F)
