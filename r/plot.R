require(tidyverse)
shapes = c(15:18, 0:4, 7:13)

otheme <- function(strip.size = 8, margin = c(.5,.5,.5,.5),
                   legend.pos = '', legend.border = F, 
                   xticks = F, yticks = F, xtitle = F, ytitle = F, 
                   xtext = F, ytext = F, xgrid = F, ygrid = F,
                   xsize = 8, ysize = 8) {
    #{{{ custom theme
    o = theme_bw() +
        theme(strip.background = element_blank(),
              strip.text = element_text(size = strip.size,
                                        margin = margin(0,0,0,0,'lines')))
    #{{{ legend
    o = o + 
        theme(legend.background = element_blank(),
              legend.title = element_blank(),
              legend.key.size = unit(.8, 'lines'), 
              legend.text = element_text(size = 8))
    if(legend.pos == 'none') {
        o = o + theme(legend.position = 'none')
    } else if(legend.pos == 'top.center.out') {
        o = o + theme(legend.position = c(.5,1), legend.justification = c(.5,0))
        margin[1] = margin[1] + 1 
    } else if(legend.pos == 'top.center') {
        o = o + theme(legend.position = c(.5,1), legend.justification = c(.5,1))
    } else if(legend.pos == 'top.left') {
        o = o + theme(legend.position = c(0,1), legend.justification = c(0,1))
    } else if(legend.pos == 'top.right') {
        o = o + theme(legend.position = c(1,1), legend.justification = c(1,1))
    } else if(legend.pos == 'bottom.right') {
        o = o + theme(legend.position = c(1,0), legend.justification = c(1,0))
    } else if(legend.pos == 'bottom.left') {
        o = o + theme(legend.position = c(0,0), legend.justification = c(0,0))
    }
    if(legend.border == T)
        o = o + theme(legend.box.background = element_rect())
    #}}}
    #{{{ axis tick/text/title
    if(xtitle) {
        o = o + theme(axis.title.x = element_text(size = 9))
    } else {
        o = o + theme(axis.title.x = element_blank())
    }
    if(ytitle) {
        o = o + theme(axis.title.y = element_text(size = 9))
    } else {
        o = o + theme(axis.title.y = element_blank())
    }
    if(xtext) {
        o = o + theme(axis.text.x = element_text(size = xsize))
    } else {
        o = o + theme(axis.text.x = element_blank())
    }
    if(ytext) {
        o = o + theme(axis.text.y = element_text(size = ysize))
    } else {
        o = o + theme(axis.text.y = element_blank())
    }
    if(!xticks)
        o = o + theme(axis.ticks.x = element_blank())
    if(!yticks)
        o = o + theme(axis.ticks.y = element_blank())
    #}}}
    #{{{ grid
    o = o + theme(panel.grid.minor = element_blank())
    if(!xgrid)
        o = o + theme(panel.grid.major.x = element_blank())
    if(!ygrid)
        o = o + theme(panel.grid.major.y = element_blank())
    #}}}
    o + theme(plot.margin = unit(margin, "lines"))
    #}}}
}
