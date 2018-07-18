require(tidyverse)
dirw = '/home/springer/zhoux379/data/genome/B73/corncyc'

fi = file.path(dirw, "query-results.tsv")
ti = read_tsv(fi, col_names = c('pathway', 'gids'))

ptn = "([&;\"])|(</?((i)|(sub)|(sup))>)"
to = ti %>% 
    mutate(pathway = str_replace_all(pathway, ptn, ""),
           gids = str_replace_all(gids, "[()\"]", "")) %>%
    mutate(gid = str_split(gids, " ")) %>%
    select(-gids) %>%
    mutate(pid = sprintf("cc%04d", 1:length(pathway))) %>%
    unnest() %>%
    select(pid, gid, pathway)

fo = file.path(dirw, '10.tsv')
write_tsv(to, fo)
