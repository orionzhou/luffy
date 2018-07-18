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

fo = file.path(dirw, "21.raw.ped.tsv")
write.table(tt[,c("Strain", "Parentage")], fo, sep = "\t", row.names = F, col.names = F, quote = F, na = "")

### run ped.dedup.py, manual inspection, then run ped.dedup.py again
fd = file.path(dirw, "22.dedup.tsv")
td = read.table(fd, header = F, sep = "\t", as.is = T, quote = "")
nrow(td)
peds = sort(unique(td$V2))
length(peds)
write.table(peds, file.path(dirw, "23.ped.tsv"), row.names = F, col.names = F, quote = F)


# Read in Pheno-tables
fi = file.path(dirw, "10.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
strains = unique(ti$Strain)

length(strains)
strainsn = unique(ti$Strain[ti$Year >= 2004 & ti$Year <= 2015])
strainsu = strainsn[!strainsn %in% tt$Strain]
length(strainsu)
fo = file.path(dirw, "21.parents.txt")
#write(unique(tt$Parentage), fo)
fo = file.path(dirw, "21.entries.tsv")
#write.table(tt, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = "")

fp = file.path(dirw, "22.tsv")
tp = read.table(fp, header = F, sep = "\t", as.is = T, quote = "", comment.char = "")
colnames(tp) = c("Parentage", "Parent1", "Parent2")
tp = cbind(tp, Pedigree = sprintf("%s / %s", tp$Parent1, tp$Parent2))
tt2 = merge(tt, tp, by = "Parentage", all.x = T)
stopifnot(nrow(tt) == nrow(tt2))

grp = dplyr::group_by(tt2, Strain)
#tt = tt2[,-4]
tt3 = dplyr::summarise(grp, 
	Year=Year[1], Pedigree=Pedigree[1], SeedSource=SeedSource[1], 
	GenComp=GenComp[1], UniqueTraits=UniqueTraits[1], Trial=Trial[1], 
	PreviousTesting=PreviousTesting[1])
ti2 = ti[ti$Year >= 2004 & ti$Year <= 2015,]
sum(! unique(ti2$Strain) %in% tt3$Strain)

tx = data.frame(
	name = tt3$Strain, program = "URT", aliases = "", acc = "", 	
	pedigree = tt3$Pedigree, generation = tt3$GenComp, species = "G.max", 
	comments = tt3$UniqueTraits, stringsAsFactors = F)
genmap = 1:9
names(genmap) = sprintf("F%d", 1:9)
tx$generation[!tx$generation %in% names(genmap)] = 'F1'
tx$generation = genmap[tx$generation]

tx2 = tx[tx$pedigree != '',c("name", "pedigree")]
tx2 = merge(tx2, tp[,2:4], by.x = "pedigree", by.y = "Pedigree")
tx3 = data.frame(LINE_NAME=tx2$name, PARENT_1=tx2$Parent1, PARENT_2=tx2$Parent2,
	contrib_1=0.5, contrib_2=0.5, selfing_1="FN", selfing_2="FN",
	pedigree=tx2$pedigree, stringsAsFactors = F)
colnames(tx3)[8] = "Text pedigree"
fp2 = file.path(dirw, "23.pedigree.tsv")
write.table(tx3, fp2, sep = "\t", row.names = F, col.names = T, quote = F, na = "")

f25 = sprintf("%s/25.strain.xlsx", dirw)
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
options(java.parameters = "-Xmx10000m")
fi = file.path(dirw, "10.tsv")
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
for (yr in c(2004:2006)) {
	tos = to[yrs == yr,]
#	tos$Value[tos$Value == ''] = "NULL"
	t31 = reshape(tos, direction = 'wide', timevar = "Phenotype", idvar = c("Trial", "Strain"))
	colnames(t31)[-c(1:2)] = sapply(strsplit(colnames(t31)[-c(1:2)], split = "[.]"), "[", 2)
	t31 = cbind.data.frame(t31[1:2], check = 0, t31[-c(1:2)], stringsAsFactors = F)
	t31 = t31[order(t31$Trial, t31$Strain, method = 'radix'),]
	#write.table(t32, f32, sep = "\t", row.names = F, col.names = T, quote = F, na='')
	f31 = sprintf("%s/31.pheno/%s.xlsx", dirw, yr)

colnames(t31)[1:3] = c("*Trial Code", "*Line Name", "*Check")
wb <- createWorkbook()
sheet1 <- createSheet(wb, "Sheet1")
rows   <- createRow(sheet1, 1:2)
cells  <- createCell(rows, colIndex=1:2)
setCellValue(cells[[1,1]], "Phenotype Results")
setCellValue(cells[[2,1]], "Crop")
setCellValue(cells[[2,2]], "soybean")
addDataFrame(t31, sheet1, col.names = T, row.names = F, startRow = 3, startColumn = 1, showNA = T, characterNA = "NULL")
saveWorkbook(wb, f31)
}


#------obsolete
grp = dplyr::group_by(tt, Strain)
tt2 = dplyr::summarise(grp, 
	Year=paste(unique(Year[Year!='']), collapse=", "), 
	Parentage=Parentage[1], nParentage=length(unique(Parentage)), 
	SeedSource=paste(unique(SeedSource[SeedSource!='']), collapse=", "), 
	GenComp=paste(unique(GenComp[GenComp!='']), collapse=", "), 
	UniqueTraits=paste(unique(UniqueTraits[UniqueTraits!='']), collapse=", "), 
	Trial=paste(unique(Trial[Trial!='']), collapse=", "), 
	PreviousTesting=paste(unique(PreviousTesting[PreviousTesting!='']), collapse=", "))
table(tt2$nParentage)
#pnt = tt$Parentage
#idxs = grep("[^([{] ?[xX] ?[^([]", pnt)

## process all_ped data
fi = file.path(dirw, "all_ped_04_15.RData")
x3 = load(fi)
tp = unique(all_ped[,c(-2,-3)])
tp[is.na(tp)] = ''
grp = dplyr::group_by(tp, Strain)
tp = dplyr::summarise(grp, Female=Female[1], Male=Male[1], 
	Synonyms=paste(unique(Synonyms[Synonyms!='']), collapse=", "), 
	Comments=paste(unique(Comments[Comments!='']), collapse=", "))

length(strains)
strainsu = strains[!strains %in% tt$Strain]
length(strainsu)
strainsu = strains[!strains %in% tt$Strain & !strains %in% ap$Strain]
length(strainsu)
head(strainsu, 50)


tm = merge(tp, tt, by='Strain', all.x=T)
tm[is.na(tm)]=''

strains = unique(ti$Strain)
length(strains)
strains_unk = strains[! strains %in% tm$Strain]
length(strains_unk)

pnt = unique(c(tm$Male, tm$Female))
pnt = pnt[pnt!='']
pnt_unk = pnt[! pnt %in% tm$Strain]
length(pnt)
length(pnt_unk)

to1 = data.frame(Strain=unique(c(strains_unk, pnt_unk)), Female='', Male='', Synonyms='', Comments='', Year='', Parentage='', SeedSource='', GenComp='', UniqueTraits='', Trial='', PreviousTesting='')
to = rbind(tm, to1)
to = cbind(id = 1:nrow(to), to)

tos = to[,1:2]; colnames(tos) = c('pid','pstrain')
tf1 = merge(to, tos, by.x='Female', by.y='pstrain', all.x=T)
colnames(tf1)[ncol(tf1)] = 'maternal_id'
tf2 = merge(tf1, tos, by.x='Male', by.y='pstrain', all.x=T)
colnames(tf2)[ncol(tf2)] = 'paternal_id'

to = tf2[,c(3,4,1,2,5:15)]
colnames(to) = c('id','name','maternal','paternal','alias_pi','alias_experiment',
	'year','parentage','seedsource','generation','unique_traits','trial',
	'previous_testing','maternal_id','paternal_id')
to = to[order(to$id),]
fo = file.path(dirw, "11.strain.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = T, na='')

tx = data.frame(name=to$name, program="URT", aliases=to$alias_pi, acc="", pedigree=to$parentage, generation=to$generation, species="Y", comments=to$unique_traits, stringsAsFactors = F)
genmap = 1:9
names(genmap) = sprintf("F%d", 1:9)
tx$generation[!tx$generation %in% names(genmap)] = 'F1'
tx$generation = genmap[tx$generation]
fo = file.path(dirw, "11.strain.t3.tsv")
write.table(tx, fo, sep = "\t", row.names = F, col.names = T, quote = F, na='')

