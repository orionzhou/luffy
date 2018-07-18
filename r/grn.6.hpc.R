require(randomForest)
source("/home/springer/zhoux379/source/genie3/GENIE3_R/genie3.R")

dirw = file.path(Sys.getenv("misc2"), "grn23", "62.genie3")

fi = file.path(dirw, "01.matrix.tsv")
fr = '/home/springer/zhoux379/data/genome/Zmays_v4/TF/11.TF.txt'
rids = scan(fr, what = character())

expr.matrix <- read.expr.matrix(fi, form="rows.are.samples")
gids = rownames(expr.matrix)
rids = rids[rids %in% gids]

weight.matrix <- get.weight.matrix(expr.matrix, input.idx=rids)
link.list <- get.link.list(weight.matrix, report.max=10000)

save(weight.matrix, link.list, file = file.path(dirw, "10.genie3.rda"))