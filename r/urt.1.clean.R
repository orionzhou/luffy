require(dplyr)
require(ggmap)
require(xlsx)
require(lubridate)

dirw = file.path(Sys.getenv("misc2"), "urt")

### read raw Liana data
yrs = 1989:2015
udr = list()
for (yr in yrs) {
	fi = sprintf("%s/raw_liana/%d_out/urt_data.RData", dirw, yr)
	x = load(fi)
	for (skey in names(urt_data)) {
		nkey = sprintf("Y%d_%s", yr, skey)
		udr[[nkey]] = urt_data[[skey]]
	}
}

### read processed (final) Liana data
fi = file.path(dirw, "urt_data.RData")
x = load(fi)
ud = urt_data; pl = planting_maturity; yp = yield_plots

### re-read Maturity tables
for (tname in names(pl)) {
	tm = udr[[tname]]
	
	ps = strsplit(tname, split="[_]")
	yr = as.integer(substr(ps[[1]][1], 2, 5))
	
	rnames = toupper(trimws(tm[,1]))
	ridx1 = which(rnames == 'DATEPLANTED')
	ridx2 = which(rnames == 'DAYSTOMATURE')
	stopifnot(length(ridx1) == 1 && length(ridx2) == 1)
	numcol = ncol(tm)
	
	flag = 1
	dp = as.Date(as.character(tm[ridx1, 2:numcol]), "%d-%b")
	if(sum(!is.na(dp)) == 0) {
		dp = as.Date(as.character(tm[ridx1, 2:numcol]), "%m/%d")
		flag = 2
	}
	stopifnot(sum(!is.na(dp)) > 0)
	year(dp) = yr
	dm = as.integer(tm[ridx2, 2:numcol])
	
	ridxl = ridx1 - 1
	if(flag == 1) {
		col2 = as.Date(as.character(tm[1:ridxl,2]), "%d-%b")
	} else {
		col2 = as.Date(as.character(tm[1:ridxl,2]), "%m/%d")
	}
	ridxf = which(!is.na(col2))
	
	if(flag == 1) {
		rowf = as.Date(as.character(tm[ridxf,2:numcol]), "%d-%b")
	} else {
		rowf = as.Date(as.character(tm[ridxf,2:numcol]), "%m/%d")
	}
	year(rowf) = yr
	
	to = tm[1:ridxl,]
	to[ridxf,2:numcol] = as.character(rowf)
	for (ridx in 1:ridxl) {
		if (ridx == ridxf) next
		to[ridx,2:numcol] = as.character(rowf + as.integer(to[ridx,]))
	}
	pl1 = data.frame(DatePlanted = as.character(dp), DaystoMature = dm, stringsAsFactors = F)
	rownames(pl1) = colnames(tm)[2:numcol]
	ud[[tname]] = to
	pl[[tname]] = pl1
}

n_nns = sapply(pl, myfunc <- function(x) sum(!is.na(x[,'DatePlanted'])))
tnames_fix = names(n_nns)[n_nns == 0]
stopifnot(length(tnames_fix) == 0)

save(ud, pl, yp, file = file.path(dirw, "00.urt.rda"))

