strjoin = function(x, sep=" ") paste(c(x), sep=sep, collapse=sep)

strconcat <- function(x) { paste(x, sep=" ", collapse=" ") }

stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, colour="red", geom=geom, width=0.2, ...)
}

rowConvert <- function(x, cs, vs) {
  y = cbind(matrix(rep(as.character(x[cs]),length(vs)), nrow=length(vs), byrow=TRUE), names(x)[vs], x[vs]);
  colnames(y) = c(names(x)[cs], 'variable', 'value')
  rownames(y) = c()
  y
}

.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           capture.output(print(object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

# shorthand
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

lg <- function(a, b, name1, name2, f_png) {
  png(filename=f_png, width=500, height=500, units='px');
  plot(a, b, type="p", xlab=name1, ylab=name2)
  fit = lm(b~a)
  abline(fit, col="blue")
  fit.sum = summary(fit)
  ann = paste("adjusted Rsquare = ", sprintf("%.04f", fit.sum$adj.r.squared), sep="")
  text(0.8*min(a)+0.2*max(a), 0.8*min(b)+0.2*max(b), ann, col='red')
  dev.off();
}

get_mt_ids <- function(opt) {
  f_id = "../conf/acc_ids.tbl"
  idt = read.table(f_id, sep="\t", header=T, stringsAsFactors=F)
  idt$id[which(idt[,opt]==1)]
}

sortC <- function(...) {
    a <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", a))
    Sys.setlocale("LC_COLLATE", "C")
    sort(...)
}

plot_hist <- function(x, fo = 'tmp.pdf', xlab='xlab', ylab='count') {
p1 = ggplot() +
  geom_histogram(aes(x = x)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.3,0.1,0.1,0.1), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
ggsave(p1, filename = fo, width = 6, height = 6)
}

