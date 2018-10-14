require(tidyverse)
require(RColorBrewer)
require(viridis)
require(ggsci)

tcd = brewer.pal.info
tc = as_tibble(tcd) %>%
    add_column(pid = rownames(tcd)) %>%
    transmute(set = category, pid = pid, numcol = maxcolors)

tp = tibble()
for (i in 1:nrow(tc)) {
    set = tc$set[i]
    numcol = tc$numcol[i]
    pid = tc$pid[i]
    pcols = brewer.pal(numcol, pid)
    tp1 = tibble(set = set, pid = pid, cid = 1:numcol, col = pcols)
    tp = rbind(tp, tp1)
}
viridis_pals = c("magma", "inferno", "plasma", "viridis", "cividis")
for (i in 1:5) {
    set = 'viridis'
    numcol = 12
    pid = viridis_pals[i]
    pcols = viridis_pal(option = pid)(numcol)
    tp1 = tibble(set = set, pid = pid, cid = 1:numcol, col = pcols)
    tp = rbind(tp, tp1)
}
sci_pals = c("npg","aaas","nejm","lancet","jama","jco","ucscgb","d3","locuszoom",
             "igv","uchicago","startrek","tron","futurama","rickandmorty",
             "simpsons",'gsea','material')
for (pid in sci_pals) {
    set = 'ggsci'
    palstr = sprintf("pal_%s()(18)", pid)
    pcols = eval(parse(text = palstr))
    pcols = pcols[!is.na(pcols)]
    numcol = length(pcols)
    tp1 = tibble(set = set, pid = pid, cid = 1:numcol, col = pcols)
    tp = rbind(tp, tp1)
}

tp = tp %>% mutate(coord = sprintf("%s_%d", pid, cid)) %>%
    mutate(coord = factor(coord, levels = coord))
p1 = ggplot(tp) +
    geom_tile(aes(x = cid, y = pid, fill = coord), color = 'white', size = 0.5) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    facet_grid(set~., scale = 'free_y', space = 'free_y') + 
    scale_fill_manual(values = tp$col) + 
    theme_bw() +
    theme(legend.position = 'none') +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8)) +
    theme(axis.text.y = element_text(size = 8, color = "black", angle = 0))
fo = file.path('/home/springer/zhoux379/git/luffy/rmd/color.palettes.pdf')
ggsave(p1, filename = fo, width = 7, height = 9)
