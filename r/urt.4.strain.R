require(dplyr)
require(ggmap)
require(xlsx)
require(plyr)

dirw = file.path(Sys.getenv("misc2"), "urt")
fi = file.path(dirw, "00.urt.rda")
x = load(fi)

## read raw entry data
fi = file.path(dirw, "URT_Entries04_05.RData")
x2 = load(fi)
et = urt_entries

tnames = names(ud)
yrs = sapply(strsplit(tnames, split="[_]"), "[", 1)
exs = sapply(strsplit(tnames, split="[_]"), "[", 2)
phs = sapply(strsplit(tnames, split="[_]"), "[", 3)
idxs = which(phs == "Entries")
tt = data.frame()
for (idx in idxs) {
	ds = ud[[idx]]; yr = yrs[idx]; test = exs[idx]
	cols = colnames(ds)
	n = nrow(ds)
	if("SeedSource" %in% cols) {seeds = ds$SeedSource} else {seeds = rep(NA, n)}
	if("PreviousTesting" %in% cols) {ptests = ds$PreviousTesting} else {ptests = rep(NA, n)}
	if("UniqueTraits" %in% cols) {utrs = ds$UniqueTraits} else {utrs = rep(NA, n)}
	tts = data.frame(Year=rep(yr, nrow(ds)), Ent=ds$Ent, Strain=ds$Strain, Parentage=ds$Parentage, SeedSource=seeds, GenComp=ds$GenComp, UniqueTraits=utrs, Trial=rep(test, nrow(ds)), PreviousTesting=ptests, stringsAsFactors=F)
	tt = rbind(tt, tts)
}

tt$Year = substr(tt$Year, 2, 5)
tt$Parentage[tt$Parentage=='na'] = ''
tt[is.na(tt)] = ''
tt$Parentage = gsub("\n", "", tt$Parentage)
tt$Parentage = gsub("\r", "", tt$Parentage)
tt$Parentage = gsub("'", "", tt$Parentage)
tt$Parentage = gsub("\"", "", tt$Parentage)
tt$Parentage = gsub("#", "", tt$Parentage)

tt2 = tt[,c('Strain', 'Year', 'Parentage', 'SeedSource', 'GenComp', 'UniqueTraits', 'Trial', 'PreviousTesting')]
colnames(tt2)[3] = 'Pedigree'
tt2 = ddply(tt, .(Strain), summarise, i=which(nchar(Parentage)==max(nchar(Parentage)))[1], 
	Year=Year[i], Pedigree=Parentage[i], SeedSource=SeedSource[i], 
	GenComp=GenComp[i], UniqueTraits=UniqueTraits[i], Trial=Trial[i], 
	PreviousTesting=PreviousTesting[i])
tt2 = tt2[,-2]

fi = file.path(dirw, "10.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
psts = unique(ti$Strain)
psts = psts[psts!='_']
nsts = psts[!psts %in% tt$Strain]
length(nsts)
tn = data.frame(Strain=nsts, Year='', Pedigree='', SeedSource='', GenComp='', UniqueTraits='', Trial='', PreviousTesting='', stringsAsFactors = F)

to = rbind(tt2, tn)
nrow(to)
#stopifnot(nrow(to) == length(unique(to$Strain)))
fo = file.path(dirw, "20.ped.raw.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = "")

fo = file.path(dirw, "21.tsv")
write.table(to[,c("Strain", "Pedigree")], fo, sep = "\t", row.names = F, col.names = F, quote = F, na = "")

## see apple notes

### update ped-file
fi = file.path(dirw, "20.ped.raw.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T, quote = "")

f_fullped = file.path(dirw, "28.dedup.tsv")
t_fullped = read.table(f_fullped, header = F, sep = "\t", as.is = T, quote = "")[,1:2]
t_fullped = t_fullped[t_fullped$V2 != '',]
colnames(t_fullped) = c('sid', 'FullPed')

fr = file.path(dirw, "27.rename.tsv")
tr = read.table(fr, header = F, sep = "\t", as.is = T, quote = "")
stopifnot(ti$Strain %in% tr$V3)
dicm = tr$V1; names(dicm) = tr$V3
t_alias = ddply(tr, .(V1), summarise, aliases = paste(V3, sep = ","))
dica = t_alias$aliases; names(dica) = t_alias$V1

fx = file.path(dirw, "29.tsv")
tx = read.table(fx, header = F, sep = "\t", as.is = T, quote = "")
#stopifnot(tr$V1 %in% tx$V1)

ti2 = ti
ti2$Strain = dicm[ti2$Strain]
ti3 = ddply(ti2, .(Strain), summarise, 
	Year=Year[1], SeedSource=SeedSource[1], 
	GenComp=GenComp[1], UniqueTraits=UniqueTraits[1], Trial=Trial[1], 
	PreviousTesting=PreviousTesting[1])
tx2 = merge(tx, ti3, by.x = 'V1', by.y = 'Strain', all.x = T)
stopifnot(nrow(tx2)==nrow(tx))

## update observation table
fv = file.path(dirw, "10.tsv")
tv = read.table(fv, header = T, sep = "\t", as.is = T)
tv = tv[tv$Strain != "_",]
stopifnot(tv$Strain %in% tr$V3)

tv$Strain = dicm[tv$Strain]
fo = file.path(dirw, "30.observations.tsv")
write.table(tv, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = "")

### make t3-style strain table
ti = tx2
colnames(ti)[1:4] = c("Strain", "Parent1", "Parent2", "Pedigree")
tx = data.frame(
	name = ti$Strain, program = "URT", aliases = "", acc = "", 	
	pedigree = ti$Pedigree, generation = ti$GenComp, species = "G.max", 
	comments = ti$UniqueTraits, stringsAsFactors = F)
genmap = 1:9
names(genmap) = sprintf("F%d", 1:9)
tx$generation[!tx$generation %in% names(genmap)] = 'F1'
tx$generation = genmap[tx$generation]
#tx$aliases = dica[tx$name]

tx3 = data.frame(LINE_NAME=ti$Strain, 
	PARENT_1=ti$Parent1, PARENT_2=ti$Parent2,
	contrib_1=0.5, contrib_2=0.5, selfing_1="FN", selfing_2="FN",
	pedigree=ti$Pedigree, stringsAsFactors = F)
colnames(tx3)[8] = "Text pedigree"
tx4 = merge(tx3, t_fullped, by.x = 'LINE_NAME', by.y = 'sid', all.x = T)
fp = file.path(dirw, "30.allstrains.ped.txt")
write.table(tx3, fp, sep = "\t", row.names = F, col.names = T, quote = F, na = "")

f25 = sprintf("%s/30.allstrains.xlsx", dirw)
wb <- createWorkbook()
sheet1 <- createSheet(wb, "Sheet1")
rows   <- createRow(sheet1, 1:5)
cells  <- createCell(rows, colIndex=1:8)
setCellValue(cells[[1,1]], "Line Submission Form")
setCellValue(cells[[2,1]], "Version")
setCellValue(cells[[2,2]], "30Sep16")
setCellValue(cells[[3,1]], "*Crop")
setCellValue(cells[[3,2]], "Soybean")
setCellValue(cells[[4,1]], "*Line Name")
setCellValue(cells[[4,2]], "*Breeding Program")
setCellValue(cells[[4,3]], "Aliases")
setCellValue(cells[[4,4]], "GRIN Accession")
setCellValue(cells[[4,5]], "Pedigree")
setCellValue(cells[[4,6]], "*Filial Generation")
setCellValue(cells[[4,7]], "*Species")
setCellValue(cells[[4,8]], "Comments")
setCellValue(cells[[5,3]], "comma separated values")
setCellValue(cells[[5,5]], "Purdy notation")
setCellValue(cells[[5,6]], "0 - 9")
setCellValue(cells[[5,7]], "G.max / G.soja")
addDataFrame(tx, sheet1, col.names = F, row.names = F, startRow = 6, startColumn = 1, showNA = F)
saveWorkbook(wb, f25)


## write phenotype table for upload
# options(java.parameters = "-Xmx10000m")
fi = file.path(dirw, "30.observations.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
tx = cbind.data.frame(ti, Trial = paste(ti$Test, ti$Year, ti$Location, sep = "_"), stringsAsFactors = F)[,c(7,4,3,6)]

fr = file.path(dirw, "11.trial.meta.tsv")
tr = read.table(fr, header = T, sep = "\t", as.is = T)
tr2 = tr[!is.na(tr$LocationMean) | !is.na(tr$C.V.), c(4,12,7:8,11)]
tr2 = cbind(tr2, ste = (tr2$L.S.D. / 1.96) ^ 2 / 2)
tr2 = tr2[,c(1,2,6,5)] 
colnames(tr2)[2:ncol(tr2)] = c("*Trial Mean", "*Std. Error", "*Replications")
tr3 = reshape(tr2, direction = "long", idvar = c("Trial"), varying = list(2:ncol(tr2)), timevar = 'Strain', times = colnames(tr2)[2:ncol(tr2)], v.names = 'Value')
tr3 = tr3[order(tr3$Trial),]
tr4 = cbind.data.frame(tr3, Phenotype = "YieldBuA", stringsAsFactors = F)[,c(1,2,4,3)]

to = rbind(tx, tr4)
to = tx
yrs = sapply(strsplit(to$Trial, split = "_"), "[", 2)
for (yr in unique(yrs)) {
	tos = to[yrs == yr,]
#	tos$Value[tos$Value == ''] = "NULL"
	t30 = reshape(tos, direction = 'wide', timevar = "Phenotype", idvar = c("Trial", "Strain"))
	colnames(t30)[-c(1:2)] = sapply(strsplit(colnames(t30)[-c(1:2)], split = "[.]"), "[", 2)
	t31 = cbind.data.frame(t30[1:2], check = 0, t30[-c(1:2)], stringsAsFactors = F)
	t31 = t31[order(t31$Trial, t31$Strain, method = 'radix'),]
	#write.table(t32, f32, sep = "\t", row.names = F, col.names = T, quote = F, na='')
	f31 = sprintf("%s/31.pheno/%s.txt", dirw, yr)
	colnames(t31)[1:3] = c("*Trial Code", "*Line Name", "*Check")
	
	cat("Phenotype Results\nCrop\tsoybean\n", file = f31)
	write.table(t31, f31, sep = "\t", row.names = F, col.names = T, quote = F, append = T, na = "")
}

# ls *.txt | sed 'p;s/\.txt/\.xls/' | xargs -n2 txt2xls.py