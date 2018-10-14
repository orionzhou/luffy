require(dplyr)
require(ggmap)
require(xlsx)

dirw = file.path(Sys.getenv("misc2"), "urt")

### extract trial meta-data
fi = file.path(dirw, "00.urt.rda")
x = load(fi)

tpl = data.frame()
for (i in 1:length(pl)) {
	zz = strsplit(names(pl)[i], split="[_]")[[1]]
	yr = zz[1]; ex = zz[2]
	t1 = pl[[i]]
	for (j in 1:ncol(t1)) { t1[,j] = as.character(t1[,j]) }
	t2 = t1[(!is.na(t1[,1]) & t1[,1] != '') | (!is.na(t1[,2]) & t1[,2] != 0),]
	t3 = cbind.data.frame(Year=yr, Test=ex, Location=rownames(t2), t2, stringsAsFactors = F)
	t4 = reshape(t3, direction='long', varying=list(4:ncol(t3)), idvar=c("Year", "Test", "Location"), timevar="Meta", v.names='Value', times=colnames(t3)[4:ncol(t3)])
	tpl = rbind(tpl, t4)
}

typ = data.frame()
for (i in 1:length(yp)) {
	zz = strsplit(names(pl)[i], split="[_]")[[1]]
	yr = zz[1]; ex = zz[2]
	t1 = yp[[i]]
	colnames(t1)[1] = "Meta"
	t2 = t1#[t1$Meta != 'LocationMean',]
	for (j in 2:ncol(t2)) { t2[,j] = as.character(t2[,j]) }
	t3 = cbind.data.frame(Year=yr, Test=ex, t2, stringsAsFactors = F)
	t4 = reshape(t3, direction='long', varying=list(4:ncol(t3)), idvar=c("Year", "Test", "Meta"), timevar="Location", v.names='Value', times=colnames(t3)[4:ncol(t3)])
	typ = rbind(typ, t4)
}
typ = typ[,c(1,2,4,3,5)]
unique(typ$Meta)

tt1 = rbind(tpl, typ)
tt2 = reshape(tt1, direction = 'wide', idvar=c("Year","Test","Location"), timevar="Meta")
colnames(tt2)[4:ncol(tt2)] = gsub("Value.", "", colnames(tt2)[4:ncol(tt2)])

locs_irri = c(
	'FayettevilleirrigatedAR',
	'PinetreeirrigatedAR'
)
tt3 = cbind.data.frame(tt2, Irrigation='no', stringsAsFactors = F)
tt3$Irrigation[tt3$Location %in% locs_irri] = 'yes'

fo = file.path(dirw, "02.trial.meta.tsv")
write.table(tt3, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

##### data check
fi = file.path(dirw, "01.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
fr = file.path(dirw, "02.trial.meta.tsv")
tr = read.table(fr, header = T, sep = "\t", as.is = T)

# check year
ti$Year = as.numeric(substr(ti$Year, 2, 5))
tr$Year = as.numeric(substr(tr$Year, 2, 5))
unique(ti$Year)
unique(tr$Year)
identical(unique(ti$Year), unique(tr$Year))

# check phenotypes
phs = unique(ti$Phenotype)
phs

ph_map = c(
		'Shattering30Sept' = 'Shattering', 
		'Shattering9Oct'   = 'Shattering',
		'Shattering3Oct'   = 'Shattering',
		'Shattering10Oct'  = 'Shattering'
)
#idxs = which(ti$Phenotype %in% names(ph_map))
#ti$Phenotype[idxs] = ph_map[ti$Phenotype[idxs]]

phs = unique(ti$Phenotype)
phs
#write(phs, file.path(dirw, "03.phenotype.txt"))

# check locations
lcs = unique(ti$Location)
lcs

lc_map = c(
	'UllinILDS' = 'UllinIL',
	'Eldorado??' = 'EldoradoIL',
	'NWBranch' = 'NorthwestBranchMD',
	'FayettevilleirrigatedAR' = 'FayettevilleAR',
	'LafayetteINLafayetteIN' = 'LafayetteIN',	
	'PinetreeirrigatedAR' = 'PineTreeAR'
)
idxs = which(ti$Location %in% names(lc_map))
ti$Location[idxs] = lc_map[ti$Location[idxs]]
lcs = unique(ti$Location)

idxs = which(tr$Location %in% names(lc_map))
tr$Location[idxs] = lc_map[tr$Location[idxs]]
lcs2 = unique(tr$Location)
lcs_more = lcs2[! lcs2 %in% lcs]
if(length(lcs_more) > 0) lcs = c(lcs, lcs_more)

pre = sapply(lcs, jj <- function(x) substr(x, start=1, stop=nchar(x)-2))
suf = sapply(lcs, jj <- function(x) substr(x, start=nchar(x)-1, stop=nchar(x)))
idxs = which(suf %in% state.abb)
ts1 = data.frame(name = lcs[idxs], city = pre[idxs], state = suf[idxs], stringsAsFactors = F)

idxs = which(! suf %in% state.abb)
lcsc = lcs[idxs]
pre = sapply(lcsc, jj <- function(x) substr(x, start=1, stop=nchar(x)-3))
suf = sapply(lcsc, jj <- function(x) substr(x, start=nchar(x)-2, stop=nchar(x)))
ts2 = data.frame(name = lcsc, city = pre, state = suf, stringsAsFactors = F)
ts2$city[ts2$name == 'Mean'] = "Mean"; ts2$state[ts2$name == 'Mean'] = "Mean"

tl2 = rbind(ts1, ts2)
tl2 = tl2[order(tl2$state, tl2$city),]

fl = file.path(dirw, "04.location.tsv")
tl = read.table(fl, sep = "\t", as.is = T, header = T)
stopifnot(tl2$name %in% tl$name)
#get_geocode(tl2, fl)
locmap = tl$location; names(locmap) = tl$name

get_geocode <- function(tl2, fo) {
	states=tl2$state; cities=tl2$city
	states[states=='QUE'] = 'QUEBEC'
	lonlat = geocode(sprintf("%s %s", cities, states))
	tl = cbind.data.frame(tl2, location = sprintf("%s_%s", tl2$city, tl2$state), longitude=lonlat[,1], latitude=lonlat[,2])
	write.table(tl, fo, sep = "\t", row.names = F, col.names = T, quote = F)
}

##
stopifnot(ti$Location %in% tl$name)
ti$Location = locmap[ti$Location]
fo = file.path(dirw, "10.tsv")
write.table(ti, fo, sep = "\t", row.names = F, col.names = T, quote = F)

##
stopifnot(tr$Location %in% tl$name)
tr$Location = locmap[tr$Location]
tr = cbind.data.frame(tr, Trial = sprintf("%s_%s_%s", tr$Test, tr$Year, tr$Location), stringsAsFactors = F)
tr = tr[,c(1:3,13,4:12)]

tit = unique(ti[,c('Year','Test','Location')])
tit = cbind.data.frame(tit, Trial = sprintf("%s_%s_%s", tit$Test, tit$Year, tit$Location), stringsAsFactors = F)
tit2 = tit[!tit$Trial %in% tr$Trial,]
tre = cbind.data.frame(tit2, DatePlanted = sprintf("%d-01-01", tit2$Year), DaystoMature = 1, C.V. = NA, L.S.D. = NA, RowSp. = NA, Rows.Plot = NA, Reps = NA, LocationMean = NA, Irrigation = 'no', stringsAsFactors = F)

tr2 = rbind(tr, tre)
tr2 = tr2[order(tr2$Year, tr2$Test, tr2$Location),]
fo = file.path(dirw, "11.trial.meta.tsv")
write.table(tr2, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

### write trial table for upload
fr = file.path(dirw, "11.trial.meta.tsv")
tr = read.table(fr, header = T, sep = "\t", as.is = T)
fl = file.path(dirw, "04.location.tsv")
tl = read.table(fl, header = T, sep = "\t", as.is = T)

tr2 = merge(tr, tl[,c('location','longitude','latitude')], by.x = 'Location', by.y = 'location')
stopifnot(nrow(tr) == nrow(tr2))

desc = sprintf("RowSp.: %d, Rows/Plot: %d", tr2$RowSp., tr2$Rows.Plot)
pdates = as.Date(tr2$DatePlanted, "%Y-%m-%d")
hdates = pdates + tr2$DaystoMature
pdates = format(pdates, "%m/%d/%Y")
hdates = format(hdates, "%m/%d/%Y")
pdates[is.na(pdates)] = sprintf('1/1/%d', tr2$Year[is.na(pdates)])
hdates[is.na(hdates)] = sprintf('1/1/%d', tr2$Year[is.na(hdates)])
pdates = sub("0?(.+)/0?(.+)/(.+)", "\\1/\\2/\\3", pdates)
hdates = sub("0?(.+)/0?(.+)/(.+)", "\\1/\\2/\\3", hdates)
tr8 = cbind.data.frame(tr2[,c(4,2:3,1,14:15)], collab = 'Aaron Lorenz', trialdesc = '', pdate = pdates, hdate = hdates, weather = '', green = 'no', rate = '', design = '', nentry = '', nrep = '', plotsize = '', harea = '', irri = tr2$Irrigation, other = desc)
tr9 = t(tr8)
colnames(tr9) = sprintf("T%04d", 1:ncol(tr9))
rownames(tr9) = c(
"Trial Name", 
"Trial Year", 
"Experiment Name", 
"Location", 
"Latitude of field", 
"Longitude of field", 
"Collaborator", 
"Trial description", 
"Planting date", 
"Harvest date", 
"Begin weather date", 
"Greenhouse trial? (yes or no)", 
"Seeding rate (seeds/m2)", 
"Experimental design", 
"Number of entries", 
"Number of replications", 
"Plot size (m2)", 
"Harvested area (m2)", 
"Irrigation (yes or no)", 
"Other remarks")
tr9[,1:2]

psize = 2200
npart = ceiling(ncol(tr9)/psize)
for (i in 1:npart) {
fo = sprintf("%s/11.trial.meta/trial.%d.xlsx", dirw, i)
to = tr9[,((i-1)*psize+1):min(i*psize, ncol(tr9))]
wb <- createWorkbook()
sheet1 <- createSheet(wb, "Sheet1")
rows   <- createRow(sheet1, 1:4)
cells  <- createCell(rows, colIndex=1:2)
setCellValue(cells[[1,1]], "Trial Submission Form")
setCellValue(cells[[2,1]], "Template version")
setCellValue(cells[[2,2]], "20Dec2016")
setCellValue(cells[[3,1]], "Crop")
setCellValue(cells[[3,2]], "soybean")
setCellValue(cells[[4,1]], "Breeding Program Code")
setCellValue(cells[[4,2]], "URT")
addDataFrame(to, sheet1, col.names = T, row.names = T, startRow = 5, startColumn = 1)
#write.xlsx(x = tr9[,1:20], file = fo, sheetName = "Sheet1", row.names = T, col.names = T)
saveWorkbook(wb, fo)
}
