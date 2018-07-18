require(plyr)
require(rtracklayer)
require(Cairo)
require(grid)
require(Rsamtools)
require(rbamtools)
source('comp.fun.R')

## data processing functions
prep_coord_mapping <- function(dcoo, seqinfo, gap_prop = 0.4, gap_len = 5000, fill_prop = 0.02, ods = NULL) {
  dfi = dcoo
  if(ncol(dcoo) == 3) { dfi = cbind(dcoo, srd = "+", stringsAsFactors = F) }
  colnames(dfi)[1:4] = c('chr', 'beg', 'end', 'srd')
  if(!is.null(ods)) {
  	dfi = dfi[dfi$chr %in% ods,]
  	dfi = dfi[order(factor(dfi$chr, levels=ods)),]
  }
  order.idx = order(dfi[,1], dfi[,2], dfi[,3])
  df = dfi[order.idx,]
  panels = rep(0, nrow(df))
  cnt = 1
  for (i in 1:nrow(df)) {
    if( panels[i] > 0 ) { next }
    panels[i] = cnt
    idP = df[i,1]; begP = df[i,2]; endP = df[i,3]
    for (j in (i+1):nrow(df)) {
      id = df[j,1]; beg = df[j,2]; end = df[j,3]
      if( j == i | j > nrow(df) | id != idP ) { break }
      if( end < begP ) {
        gap = begP - end - 1
      } else if( endP < beg ) {
        gap = beg - endP + 1
      } else {
        gap = 0
      }
      panelBeg = min(beg, begP)
      panelEnd = max(end, endP)
      panelLen = panelEnd - panelBeg + 1
      if( gap/panelLen <= gap_prop | gap <= gap_len ) {
        panels[j] = cnt
        idP = id; begP = panelBeg; endP = panelEnd
      } else {
        break
      }
    }
    cnt = cnt + 1
  }
  
  dcon = cbind(dfi, panel_old = panels[order(order.idx)])
  panel_unique = unique(dcon$panel_old)
  panel_mapping = data.frame(panel_old = panel_unique, 
    pan = 1:length(panel_unique))
  dcon = cbind(dcon, idx = 1:nrow(dcon))
  dcon = merge(dcon, panel_mapping, by = 'panel_old')
  dcon = dcon[order(dcon$idx),]
  dcon = dcon[, !colnames(dcon) %in% c('idx','panel_old')]
  
  tmp1 = ddply(dcon, .(pan, srd), summarise, len = sum(end - beg));
  tmp2 = ddply(tmp1, .(pan), summarise, srd = srd[which(len == max(len))][1])
  dpan = ddply(dcon, .(pan), summarise, chr = unique(chr), beg = min(beg), 
    end = max(end))
  dpan = merge(dpan, tmp2, by = 'pan')

  dlen = data.frame(id = as.character(seqnames(seqinfo)), 
    size = seqlengths(seqinfo), stringsAsFactors = F)
  dpan = merge(dpan, dlen, by.x = 'chr', by.y = 'id')
  dpan = dpan[order(dpan$pan),]
  begr = as.integer( dpan$beg - 0.02 * (dpan$end - dpan$beg) )
  endr = as.integer( dpan$end + 0.02 * (dpan$end - dpan$beg) )
  dpan = cbind(dpan, begr = begr, endr = endr)
  dpan$begr = apply(dpan, 1, function(x) max(1, as.numeric(x['begr'])))
  dpan$endr = apply(dpan, 1, function(x) 
    min(as.numeric(x['size']), as.numeric(x['endr'])))
  dpan = cbind(dpan, len = dpan$endr - dpan$begr + 1)

  len_total = sum(dpan$len)
  itv = as.integer(fill_prop * len_total)
  dpan = cbind(dpan, beg.a = itv + 1)
  if(nrow(dpan) > 1) {
    for (i in 2:nrow(dpan)) {
      dpan$beg.a[i] = dpan$beg.a[i - 1] + dpan$len[i - 1] + itv + 1
    }
  }
  
  dmap = dpan[,c('pan','chr','begr','endr','srd','len','beg.a')]
  colnames(dmap)[3:4] = c('beg','end')
  cbind(dmap, end.a = dmap$beg.a + dmap$len - 1)
}
prep_ticks <- function(dmap, tick_itv) {
  dtik = data.frame(pan = c(), pos = c())  
  for (i in 1:nrow(dmap)) {
    tick_beg = pretty(c(dmap$beg[i], dmap$end[i]))[1]
    ticks = seq(tick_beg, dmap$end[i], by = tick_itv)
    ticks = ticks[ticks >= dmap$beg[i] & ticks <= dmap$end[i]]
    if(length(ticks) == 0) { ticks = tick_beg }
    dtik = rbind(dtik, data.frame(pan = rep(dmap$pan[i], length(ticks)), 
      pos = ticks))
  }
  
  dtik = merge(dtik, dmap, by = 'pan')
  dtik = cbind(dtik, pos.a = 0)
  for (i in 1:nrow(dtik)) {
    if(dtik$srd[i] == "-") {
      dtik$pos.a[i] = dtik$end.a[i] - (dtik$pos[i] - dtik$beg[i])
    } else {
      dtik$pos.a[i] = dtik$beg.a[i] + (dtik$pos[i] - dtik$beg[i])
    }
  }
  dtik = dtik[dtik$pos.a > 0, c('pan','pos','srd','pos.a')]
  dtik
}

coord_mapping <- function(dcoo, dmap) {
  if(is.null(dcoo)) { return(NULL) }
  if(nrow(dcoo)==0) { return(NULL) }
  if(ncol(dcoo) == 3) { dcoo = cbind(dcoo, srd = "+", stringsAsFactors = F) }
  colnames(dcoo)[1:4] = c('chr','beg','end','srd')
  idxs.raw = c()
  pan.idxs = c()
  for (i in 1:nrow(dmap)) {
    chr = dmap$chr[i]; beg = dmap$beg[i]; end = dmap$end[i]; pan = dmap$pan[i]
    idxss = which( dcoo$chr == chr & ( (beg <= dcoo$beg & dcoo$beg <= end) | 
      (beg <= dcoo$end & dcoo$end <= end) | 
      (beg > dcoo$beg & end < dcoo$end) ) )
    if(length(idxss) > 0) {
      idxs.raw = c(idxs.raw, idxss)
      pan.idxs = c(pan.idxs, rep(i, length(idxss)))
    }
  }
  
  if(length(idxs.raw) == 0) {idxs=NULL} else {idxs=sort(unique(idxs.raw))}
  pans = c(); begs.a = c(); ends.a = c(); srds.a = c()
  for (i in idxs) {
    chr = dcoo[i,1]; beg = dcoo[i,2]; end = dcoo[i,3]; srd = dcoo[i,4]
    
    pan.idxss = pan.idxs[ which(idxs.raw == i) ]
    pan.idx = pan.idxss[1]
    if( length(pan.idxss) > 1 ) {
      lens_ovlp = c()
      for (i in 1:length(pan.idxss) ) {
        len_ovlp = min(dmap$end[pan.idxss[i]], end) - 
          max(dmap$beg[pan.idxss[i]], beg) + 1
        lens_ovlp = c(lens_ovlp, len_ovlp)
      }
      pan.idx = pan.idxss[ which(lens_ovlp == max(lens_ovlp))[1] ]
    }
    
    pan = dmap$pan[pan.idx]
    pan.beg = dmap$beg[pan.idx]
    pan.end = dmap$end[pan.idx]
    pan.srd = dmap$srd[pan.idx]
    pan.beg.a = dmap$beg.a[pan.idx]
    pan.end.a = dmap$end.a[pan.idx]
    
    beg = max(beg, pan.beg)
    end = min(pan.end, end)
    beg.a = ifelse( pan.srd == "-", pan.end.a - (end - pan.beg), 
      pan.beg.a + (beg - pan.beg) );
    end.a = ifelse( pan.srd == "-", pan.end.a - (beg - pan.beg), 
      pan.beg.a + (end - pan.beg) );
    srd.a = ifelse(pan.srd == srd, "+", "-")
    
    pans = c(pans, pan)
    begs.a = c(begs.a, beg.a)
    ends.a = c(ends.a, end.a)
    srds.a = c(srds.a, srd.a)
  }
  cbind(dcoo[idxs,], pan = pans, beg.a = begs.a, end.a = ends.a, 
    srd.a = srds.a, stringsAsFactors = F)
}

coord_mapping_bw <- function(fbw, dmap) {
  dm = data.frame()
  
  pans = c(); begs.a = c(); ends.a = c(); srds.a = c()
  for (i in 1:nrow(dmap)) {
    chr = dmap$chr[i]; beg = dmap$beg[i]; end = dmap$end[i]; srd = dmap$srd[i]
    grs = GRanges(seqnames = chr, ranges = IRanges(beg, end = end))
    
    pan = dmap$pan[i]
    pan.beg = dmap$beg[i]; pan.end = dmap$end[i]; pan.srd = dmap$srd[i]
    pan.beg.a = dmap$beg.a[i]; pan.end.a = dmap$end.a[i]
    
    vals = import.bw(fbw, which = grs, as = "GRanges")
    begs = start(vals); ends = end(vals); scores = score(vals)
    if(length(vals) > 0) {
      if(pan.srd == "-") {
        begs.a = pan.end.a - (ends - pan.beg)
        ends.a = pan.end.a - (begs - pan.beg)
      } else {
        begs.a = pan.beg.a + (begs - pan.beg)
        ends.a = pan.beg.a + (ends - pan.beg)
      }
      srd.a = pan.srd
      dms = data.frame(chr = chr, beg = begs, end = ends, val = scores, 
        pan = pan, beg.a = begs.a, end.a = ends.a, srd.a = pan.srd, 
        stringsAsFactors = F)
      dm = rbind(dm, dms)
    }
  }
  dm
}

prep_plot_data <- function(gro, cfgs, tname, qnames, tracks, largescale = F, fill_prop = 0.02, qods = NULL) {
  tcfg = cfgs[[tname]]
  tmap = prep_coord_mapping(granges2df(gro), tcfg$seqinfo, fill_prop = fill_prop)
  gr = GRanges(seqnames = tmap$chr, ranges = IRanges(tmap$beg, end = tmap$end),
    seqinfo = tcfg$seqinfo)
    
  pres = list()
  max_len = tmap$end.a + tmap$beg.a[1] - 1
  max_pan_len = max(tmap$len)
  for (qname in qnames) {
    cfg = cfgs[[qname]]
    if(largescale) {
    	aln = read_gal(cfg$tgal, gr)
    } else {
    	aln = read_gax(cfg$tgax, gr)
    }
    if(is.null(aln)) {
      dats[[qname]] = NULL
      next
    }
    if(largescale) {
    	qmap = prep_coord_mapping(aln$tc[,7:10], cfg$seqinfo, ods=qods[[qname]], fill_prop = fill_prop)
    } else {
    	qmap = prep_coord_mapping(aln$tg[,5:8], cfg$seqinfo, ods=qods[[qname]], fill_prop = fill_prop)
    }
    grq = GRanges(seqnames = qmap$chr, 
      ranges = IRanges(qmap$beg, end = qmap$end), seqinfo = cfg$seqinfo)
    max_len = max(max_len, qmap$end.a + qmap$beg.a[1] - 1)
    max_pan_len = max(max_pan_len, qmap$Len)
    
    if(largescale) {
    	to = aln$tc
    } else {
    	to = aln$tg
    }
    pres[[qname]] = list(to = to, qmap = qmap, gr = grq)
  }
  
  tick_itv = diff( pretty(c(1, max_pan_len))[1:2] )      
      
  ttik = prep_ticks(tmap, tick_itv)
  tgap = coord_mapping(read_tabix(tcfg$gapz, gr), tmap)
  if("tgene" %in% tracks) {
  	tgene = read_tabix(tcfg$genez, gr)
  } else {
  	tgene = NULL
  }
  if(!is.null(tgene)) {
    colnames(tgene) = c('chr', 'beg', 'end', 'srd', 'id', 'type', 'cat')
    tgene = coord_mapping(tgene, tmap)
  }
  if('tmapp' %in% tracks) tmapp = coord_mapping_bw(tcfg$mapp, tmap)
  if('trnaseq' %in% tracks) {
    tr = read_bam(tcfg$rnaseq, gr, pileup = T)
    if(nrow(tr) > 0) {
      tr2 = coord_mapping(tr[,c('chr','beg','end','srd')], tmap)
      trnaseq = merge(tr, tr2, by = c('chr', 'beg', 'end', 'srd'))
    } else {
      trnaseq = NULL
    }
  }

  dats = list()
  for (qname in qnames) {
    cfg = cfgs[[qname]]
    pre = pres[[qname]]
    if(is.null(pre)) next
    to = pre$to; qmap = pre$qmap; grq = pre$gr
    
    qtik = prep_ticks(qmap, tick_itv)
    qgap = read_tabix(cfg$gapz, grq)
    qgap = coord_mapping(qgap, qmap)
    qgene = NULL
    if('qgene' %in% tracks) {
    qgene = read_tabix(cfg$genez, grq)
    if(!is.null(qgene)) {
      colnames(qgene) = c('chr', 'beg', 'end', 'srd', 'id', 'type', 'cat')
      qgene = coord_mapping(qgene, qmap)
    }
    }

    tdcoo = coord_mapping(to[,c('tid','tbeg','tend','tsrd')], tmap)
    qdcoo = coord_mapping(to[,c('qid','qbeg','qend','qsrd')], qmap)
    stopifnot(rownames(tdcoo) == rownames(to), rownames(qdcoo) == rownames(to))
    comp = cbind(to, 
      tbeg.a = tdcoo$beg.a, tend.a = tdcoo$end.a, tsrd.a = tdcoo$srd.a, 
      qbeg.a = qdcoo$beg.a, qend.a = qdcoo$end.a, qsrd.a = qdcoo$srd.a, 
      stringsAsFactors = F)
    
    snp = NULL
    if(max_len < 100000) {
      snp = read_tabix(cfg$tsnp, gr)
      snp = snp[snp$V8 == 1,]
      if(typeof(snp) == "list") if(nrow(snp) == 0) snp = NULL
      if(!is.null(snp)) {
    colnames(snp) = c('tid', 'tpos', 'ref', 'alt', 'qid', 'qpos', 'cid', 'lev')
    tsnp = coord_mapping(snp[,c(1,2,2)], tmap)
    qsnp = coord_mapping(snp[,c(5,6,6)], qmap)
    tsnp = tsnp[,c(1:2,6,8)]
    colnames(tsnp) = c("tid", 'tpos', 'tpos.a', 'tsrd.a')
    qsnp = qsnp[,c(1:2,6,8)]
    colnames(qsnp) = c("qid", 'qpos', 'qpos.a', 'qsrd.a')
    snp = merge(snp[,c(1:2,5:8)], tsnp, by = c('tid', 'tpos'))
    snp = merge(snp, qsnp, by = c('qid', 'qpos'))
      }
    }
    
    dat = list(tmap = tmap, grt = gr, qmap = qmap, grq = grq, 
      ttik = ttik, tgap = tgap, tgene = tgene, 
      qtik = qtik, qgap = qgap, qgene = qgene,
      comp = comp, snp = snp)
    
    if('tmapp' %in% tracks) dat[['tmapp']] = tmapp
    if('trnaseq' %in% tracks) dat[['trnaseq']] = trnaseq
    if('msnp' %in% tracks) {
      msnp = read_tabix(cfg$vsnp, gr)
      if(!is.null(msnp)) {
      colnames(msnp) = c('chr', 'pos', 'alt', 'gt', 'rd', 'qual', 'mapqual')
      msnp = msnp[msnp$gt == 2,]
      msnp = coord_mapping(msnp[,c(1,2,2)], tmap)
      msnp = cbind(msnp[,c(1:2)], pos.a = msnp$beg.a, stringsAsFactors = F)
      }
      dat[['msnp']] = msnp
    }
    if('tsnp' %in% tracks) {
      tsnp = NULL
      if(!is.null(snp)) {
      tsnp = snp[,c('tid','tpos','tpos.a')]
      colnames(tsnp) = c('chr','pos','pos.a')
      }
      dat[['tsnp']] = tsnp
    }
    if('tpct' %in% tracks) {
      tpct = coord_mapping_bw(cfg$tpct, tmap)
      dat[['tpct']] = tpct
    }
    if('mcov' %in% tracks) {
      mcov = coord_mapping_bw(cfg$vcov, tmap)
      dat[['mcov']] = mcov
    }
    if('qpacbio' %in% tracks) {
      dat[['qpacbio']] = NULL
      if(!is.null(cfg$pacbio)) {
      tr = read_bam(cfg$pacbio, grq, pileup = T)
      if(nrow(tr) > 0) {
        tr2 = coord_mapping(tr[,c('chr','beg','end','srd')], qmap)
        dat[['qpacbio']] = merge(tr, tr2, by = c('chr', 'beg', 'end', 'srd'))
      }
      }
    }
    if('qrnaseq' %in% tracks) {
      tr = read_bam(cfg$rnaseq, grq, pileup = T)
      dat[['qrnaseq']] = NULL
      if(nrow(tr) > 0) {
        tr2 = coord_mapping(tr[,c('chr','beg','end','srd')], qmap)
        dat[['qrnaseq']] = merge(tr, tr2, by = c('chr', 'beg', 'end', 'srd'))
      }
    }
    dats[[qname]] = dat
  }
  dats$max_len = max_len
  dats
}
comp.plot <- function(dats, tname, qnames, tracks, scale.ht = unit(0.8, 'npc'), draw.legend.gene = F, legend.opt = 'all', largescale = F) {
  ## some default plotting parameters
  tracktypes = c('tgene' = 'gene', 'taxis' = 'axis',  'tgap' = 'gap', 
    'qgap' = 'gap', 'qaxis' = 'axis', 'qgene' = 'gene',
    'link' = 'link', 'tmapp' = 'mapp', 
    'mcov' = 'mcov', 'msnp' = 'snp', 'tsnp' = 'snp', 'tpct' = 'tpct', 
    'qpacbio' = 'bam', 'qrnaseq' = 'bam', 'trnaseq' = 'bam')

  trackheight = c('axis' = 30, 'gap' = 10, 'gene' = 15, 'link' = 45,
    'snp' = 10, 'mcov' = 20, 'mapp' = 20, 'tpct' = 20, 'bam' = 30)
  if(largescale) {trackheight['link'] = 90}

  fillg = c('TE' = 'slategray3', 'Gene' = 'brown4', 
    'CC-NBS-LRR' = 'forestgreen', 'TIR-NBS-LRR' = 'forestgreen', 
    'CRP0000-1030' = 'dodgerblue', 'NCR' = 'dodgerblue', 'CRP1600-6250' = 'dodgerblue')

  left.wd = 80 # left (legend) panel width
  top.ht = 15 # top (title) panel height (if to draw)
  top.ht = ifelse(draw.legend.gene, top.ht, 0)

  max_len = dats$max_len
  main = sprintf("%s compare to %d accessions", tname, length(qnames))
  qheight = sum(trackheight[tracktypes[tracks]])
  height = top.ht + length(qnames) * qheight
  
  idx_cmp = which(tracks == 'link')
  qhtb = 0; qhte = sum(trackheight[tracktypes[tracks[1:(idx_cmp-1)]]])
  thtb = sum(trackheight[tracktypes[tracks[1:idx_cmp]]]); thte = qheight
  
  vtr <- viewport(x = unit(left.wd, 'points'), y = unit(1, 'npc'), 
    width = unit(1, 'npc') - unit(left.wd, 'points'), 
    height = unit(top.ht, 'points'), xscale = c(1, max_len), 
    just = c('left', 'top'), name = 'headr')
  
  vbl <- viewport(x = 0, y = unit(1, 'npc') - unit(top.ht, 'points'), 
    width = unit(left.wd, 'points'), 
    height = unit(1, 'npc') - unit(top.ht, 'points'), 
    just = c('left', 'top'), name = 'left')

  vbr <- viewport(x = unit(left.wd, 'points'), 
    y = unit(1, 'npc') - unit(top.ht, 'points'), 
    width = unit(1, 'npc') - unit(left.wd, 'points'), 
    height = unit(1, 'npc') - unit(top.ht, 'points'), 
    xscale = c(1, max_len), 
    just = c('left', 'top'), name = 'right')
  
  vl <- viewport(x = 0, y = unit(1, 'npc'), 
    width = unit(left.wd, 'points'), height = unit(1, 'npc'), 
    just = c('left', 'top'), name = 'left')
  
  vr <- viewport(x = unit(left.wd, 'points'), y = unit(1, 'npc'), 
    width = unit(1, 'npc') - unit(left.wd, 'points'), height = unit(1, 'npc'), 
    xscale = c(1, max_len), 
    just = c('left', 'top'), name = 'right')
  
  if(draw.legend.gene) {
    grobs = gList(plot_legend_gene(x = unit(0, 'npc'), opt = legend.opt, vp = vtr))
#    grob_title = plot_title(main = '', subtitle = '', vp = vtr)
  } else {grobs = gList()}
  
  os = top.ht
  for (name in qnames) {
    pht = sum(trackheight[tracktypes[tracks]])
    grob_q = rectGrob(x = 0, width = 1, 
      y = unit(1, 'npc') - unit(os + qhtb, 'points'), 
      height = unit(qhte-qhtb, 'points'), just = c('left', 'top'),
      gp = gpar(fill = 'lemonchiffon', alpha = 1, col = NA, lwd = 1), vp = vl)
    grob_t = rectGrob(x = 0, width = 1, 
      y = unit(1, 'npc') - unit(os + thtb, 'points'), 
      height = unit(thte-thtb, 'points'), just = c('left', 'top'),
      gp = gpar(fill = 'lavenderblush', alpha = 1, col = NA, lwd = 1), vp = vl)
    grob_panel = rectGrob(x = 0, 
      y = unit(1, 'npc') - unit(os, 'points'), 
      width = 1, height = unit(pht, 'points'), just = c('left', 'top'),
      gp = gpar(fill = NA, alpha = 1, col = 'black', lwd = 0.5))
    grobs = gList(grobs, grob_q, grob_t, grob_panel)
    os = os + pht
  }
  grob_grid = plot_grid(vp = vbr)
  grobs = gList(grobs, grob_grid)
  
  tgt.up = which(tracks == 'taxis') < which(tracks == 'link')
  os = top.ht
  for (name in qnames) {
    dat = dats[[name]]
    pht = sum(trackheight[tracktypes[tracks]])
        
    oos = 0
    for (track in tracks) {
      tracktype = tracktypes[track]
      tht = trackheight[tracktype]
      y2 = unit(1, 'npc') - unit(os + oos, 'points')
      y1 = unit(1, 'npc') - unit(os + oos + tht, 'points')
      y = unit(1, 'npc') - unit(os + oos + tht/2, 'points')
      if(tracktype == 'axis') {
        if(track == 'taxis') {dmap = dat$tmap} else {dmap = dat$qmap}
        if(track == 'taxis') {dtik = dat$ttik} else {dtik = dat$qtik}
        lab = ifelse(track == 'taxis', tname, name)
        grob_le = plot_legend(sprintf("%s", lab), y = y, vp = vl)
        grob_ri = plot_axis(dmap, dtik, y = y, text.above = T, vp = vr, largescale = largescale)
      } else if(tracktype == "gap") {
        if(track == 'tgap') {dgap = dat$tgap} else {dgap = dat$qgap}
        lab = ifelse(track == 'tgap', tname, name)
        grob_le = plot_legend(sprintf("%s gap", lab), y = y, vp = vl)
        grob_ri = plot_gap(dgap, y = y, vp = vr)
      } else if(tracktype == "gene") {
        if(track == 'tgene') {dgen = dat$tgene} else {dgen = dat$qgene}
        lab = ifelse(track == 'tgene', tname, name)
        grob_le = plot_legend(sprintf("%s gene", lab), y = y, vp = vl)
        grob_ri = plot_gene(dgen, fillg, y = y, vp = vr)
      } else if(tracktype == "link") {
        grob_le = plot_legend(sprintf("%s/%s", name, tname), y = y, vp = vl)
        grob_ri = plot_link(dat$comp, dat$snp, y = y, height = tht, tgt.up = tgt.up, vp = vr)
      } else if(tracktype == "snp") {
        stopifnot(track %in% c("msnp", "tsnp"))
        lab = ifelse(track == 'msnp', "Mapping", "Synteny")
        if(track == 'msnp') {snp = dat$msnp} else {snp = dat$tsnp}
        grob_le = plot_legend(sprintf("%s SNP", lab), y = y, vp = vl)
        grob_ri = plot_snp(snp, y = y, vp = vr)
      } else if(tracktype == "bam") {
        stopifnot(track %in% c("qpacbio", "qrnaseq", "trnaseq"))
        if(track == "qpacbio") {
          lab = "PacBio"; tr = dat$qpacbio
        } else if(track == 'qrnaseq') {
          lab = "RNA-Seq"; tr = dat$qrnaseq
        } else {
          lab = "RNA-Seq"; tr = dat$trnaseq
        }
        grob_le = plot_legend(lab, y = y, vp = vl)
        if(is.null(tr)) { grob_ri = NA; next }
        tr = tr[tr$y <= 20,]
        vd <- viewport(x = unit(left.wd, 'points'), y = y, 
          width = unit(1, 'npc') - unit(left.wd, 'points'), 
          height = unit(tht - 3, 'points'), 
          xscale = c(1, max_len), 
          yscale = c(min(tr$y), max(tr$y)),
          just = c('left', 'center'))
        grob_ri = plot_reads(tr, vp = vd)
      } else if(tracktype == "mcov") {
        grob_le = plot_legend("MappingCoverage", y = y, vp = vl)
        vd <- viewport(x = unit(left.wd, 'points'), y = y, 
          width = unit(1, 'npc') - unit(left.wd, 'points'), 
          height = unit(tht - 3, 'points'), 
          xscale = c(1, max_len), 
          yscale = c(min(dat$mcov$val), max(dat$mcov$val)),
          just = c('left', 'center'))
        grob_ri = plot_hist(dat$mcov, vp = vd)
      } else if(tracktype == "tpct") {
        grob_le = plot_legend("SNPdensity", y = y, vp = vl)
        vd <- viewport(x = unit(left.wd, 'points'), y = y, 
          width = unit(1, 'npc') - unit(left.wd, 'points'), 
          height = unit(tht - 3, 'points'), 
          xscale = c(1, max_len), 
          yscale = c(0, 1),
          just = c('left', 'center'))
        grob_ri = plot_hist_rect(dat$tpct, vp = vd)
      } else if(tracktype == "mapp") {
        lab = ifelse(track == 'tmapp', tname, name)
        grob_le = plot_legend("Mappability", y = y, vp = vl)
        stopifnot(track == 'tmapp')
        mapp = dat$tmapp
        vd <- viewport(x = unit(left.wd, 'points'), y = y, 
          width = unit(1, 'npc') - unit(left.wd, 'points'), 
          height = unit(tht - 3, 'points'), 
          xscale = c(1, max_len), yscale = c(0, 1),
          just = c('left', 'center'))
        grob_ri = plot_hist(dat$tmapp, vp = vd)
      } else {
        stop(cat("unknonwn tracktype", tracktype, track, "\n"))
      }
      grobs = gList(grobs, grob_le, grob_ri)
      oos = oos + tht
    }
    os = os + pht
  }
  grobs = gList(grobs, plot_scale(max_len, y = scale.ht, vp = vbr, largescale = largescale))
  list(ht = height, grobs = grobs)
}

## plotting functions
plot_grid <- function(wd = 1000, itv = 30, vp = NULL) {
  n = wd %/% itv
  cols = rep('lightsteelblue3', n + 1)
  cols[1] = 'red'
  xs = unit(seq(0, length.out = n + 1, by = itv), 'points')
  segmentsGrob( x0 = xs, x1 = xs, y0 = 0, y1 = 1, 
    gp = gpar(col = cols, alpha = 1, lwd = 0.5), vp = vp)
}
plot_legend <- function(txt, x = unit(4, 'points'), y = unit(0.5, 'npc'), 
  text.rot = 0, vp = NULL) {
  txtmap = c("HM101" = "A17 (Mt4.0)", "HM340.FN" = "R108 (v1.0)", "HM340.FN/HM101" = "", "HM340" = "R108 (ALLPATHS)", "HM340.PB" = "R108 (PacBio)", "HM340/HM101" = "", "HM340.PB/HM101" = "")
  if(txt %in% names(txtmap)) {txt = txtmap[txt]}
  textGrob( label = txt, 
    x = x, y = y, just = c("left", "center"), 
    rot = text.rot, gp = gpar(cex = 0.7, fontfamily = "serif"), 
    vp = vp)
}
plot_axis <- function(dmap, dtik, y = unit(0.5, 'npc'), 
  col.p = 'red', col.n = 'blue', text.above = F, text.rot = 0, largescale = F, vp = NULL) {
  if(is.null(dmap)) return(NULL)
  dax = dmap[,c('chr', 'beg.a', 'end.a', 'srd')]
  colnames(dax) = c('id','beg','end','srd')
  
  line.cols = rep(col.p, nrow(dax))
  line.cols[which(dax$srd == "-")] = col.n
  axisgrob = segmentsGrob(
    x0 = unit(dax$beg, 'native'), x1 = unit(dax$end, 'native'),
    y0 = y, y1 = y, 
    gp = gpar(col = line.cols), vp = vp)
  
  text.x = c()
  text.just = c()
  for (i in 1:nrow(dax)) {
    beg = dax[i,2]; end = dax[i,3]; srd = dax[i,4]
    if(srd == '-') {
      x0 = unit(beg, 'native') + unit(5,'points')
      x1 = unit(beg, 'native')      
      y0 = y + unit(ifelse(text.above, -3, 3), 'points')
      y1 = y
      text.x = c(text.x, end)
      text.just = c(text.just, 'right', 'center')
    } else {
      x0 = unit(end, 'native') - unit(5,'points')
      x1 = unit(end, 'native')
      y0 = y + unit(ifelse(text.above, -3, 3), 'points')
      y1 = y
      text.x = c(text.x, beg)
      text.just = c(text.just, 'left', 'center')
    }
    if(i == 1) {
      arrows.x0 = x0; arrows.x1 = x1
      arrows.y0 = y0; arrows.y1 = y1
    } else {
      arrows.x0 = unit.c(arrows.x0, x0); arrows.x1 = unit.c(arrows.x1, x1)
      arrows.y0 = unit.c(arrows.y0, y0); arrows.y1 = unit.c(arrows.y1, y1)
    }
  }
  arrowgrob = segmentsGrob( 
    x0 = arrows.x0, x1 = arrows.x1, y0 = arrows.y0, y1 = arrows.y1,
    gp = gpar(col = line.cols), vp = vp)
  
  text.y = y + unit(ifelse(text.above, 5, -5), 'points') 
  text.offset = y + unit(ifelse(text.above, -10, 10), "points")
  chrgrob = textGrob( label = dax$id, 
    x = unit((dax$beg + dax$end) / 2, 'native'), 
    y = text.y, just = c("center", "center"), 
    rot = text.rot, gp = gpar(cex = 0.7),# fontfamily = "Mono"), 
    vp = vp)
    
  if( !text.above ) {
    tick.y = y + unit(3, 'points')
    text.y = y + unit(5, 'points')
    text.just = c('center', 'bottom')
  } else {
    tick.y = y - unit(3, 'points')
    text.y = y - unit(5, 'points')
    text.just = c('center', 'top')
  }
  if(is.null(dtik) | nrow(dtik) == 0) {
    gList(axisgrob, arrowgrob, chrgrob)
  } else {
  	if(largescale) {lab = dtik$pos / 1000000} else {lab = dtik$pos / 1000}
    gList(axisgrob, arrowgrob, chrgrob, 
  tikgrob = segmentsGrob(
    x0 = unit(dtik$pos.a, 'native'),
    x1 = unit(dtik$pos.a, 'native'),
    y0 = y, y1 = tick.y, 
    vp = vp),
  posgrob = textGrob( 
    label = lab,
    x = unit(dtik$pos.a, 'native'), 
    y = text.y,  just = text.just, 
    gp = gpar(cex = 0.6), vp = vp)
    )
  }
}
plot_gene <- function(dgen, fillg, y = unit(0.5, 'npc'), text.show = F, 
  text.offset = unit(10, 'points'), text.above = F, 
  text.rot = 0, vp = NULL) {
  if(is.null(dgen)) { return(NULL) }
  dgen = dgen[,c('beg.a', 'end.a', 'srd.a', 'type', 'cat')]
  colnames(dgen) = c('beg', 'end', 'srd', 'type', 'cat')
  dgen$cat[!dgen$cat %in% names(fillg)] = "Gene"
  
  linegrob1 = segmentsGrob(x0 = unit(min(dgen$beg), 'native'), 
    x1 = unit(max(dgen$end), 'native'), y0 = y, y1 = y, 
    gp = gpar(col = 'grey', lty = 3), vp = vp)
  
  dfp = dgen[dgen$srd == "+",]
  dfn = dgen[dgen$srd == "-",]
  height = 5
  yp = y + unit(height / 2, 'points')
  yn = y - unit(height / 2, 'points')
  
  grobs = gList(linegrob1)
  if( !empty(dfp) ) {
    for (cat in names(fillg)) {
      dfm = dfp[dfp$type == 'mrna' & dfp$cat == cat,]
      dfc = dfp[dfp$type == 'cds' & dfp$cat == cat,]
      if(!empty(dfm)) {
        cols = fillg[dfm$cat]; cols[is.na(cols)] = 'darkorchid'
        rnagrob = segmentsGrob(
          x0 = unit(dfm$beg, 'native'), x1 = unit(dfm$end, 'native'),
          y0 = yp, y1 = yp, gp = gpar(col = cols), vp = vp)
        grobs = gList(rnagrob, grobs)
      }
      if(!empty(dfc)) {
        cols = fillg[dfc$cat]; cols[is.na(cols)] = 'darkorchid'
        cdsgrob = rectGrob( 
          x = unit(dfc$beg, 'native'), y = yp,
          width = unit(dfc$end - dfc$beg, 'native'), 
          height = unit(height, 'points'), just = c('left', 'center'),
          gp = gpar(lwd = 0, fill = cols, alpha = 0.9), vp = vp)
        grobs = gList(cdsgrob, grobs)
      }
    }
  }
  
  if( !empty(dfn) ) {
    for (cat in names(fillg)) {
      dfm = dfn[dfn$type == 'mrna' & dfn$cat == cat,]
      dfc = dfn[dfn$type == 'cds' & dfn$cat == cat,]
      if(!empty(dfm)) {
        cols = fillg[dfm$cat]; cols[is.na(cols)] = 'darkorchid'
        rnagrob = segmentsGrob(
          x0 = unit(dfm$beg, 'native'), x1 = unit(dfm$end, 'native'),
          y0 = yn, y1 = yn, gp = gpar(col = cols), vp = vp)
        grobs = gList(rnagrob, grobs)
      }
      if(!empty(dfc)) {
        cols = fillg[dfc$cat]; cols[is.na(cols)] = 'darkorchid'
        cdsgrob = rectGrob( 
          x = unit(dfc$beg, 'native'), y = yn,
          width = unit(dfc$end - dfc$beg, 'native'), 
          height = unit(height, 'points'), just = c('left', 'center'),
          gp = gpar(lwd = 0, fill = cols, alpha = 0.9), vp = vp)
        grobs = gList(cdsgrob, grobs)
      }
    }
  }
  grobs
}
plot_gap <- function(dgap, y = unit(0.5, 'npc'), fill = 'grey', vp = NULL) {
  if(is.null(dgap)) { return(NULL) }
  dgap = dgap[,c('chr', 'beg.a', 'end.a')]
  colnames(dgap) = c('id','beg','end')
  ht = unit(5, 'points')
  rectGrob(
    x = unit(dgap$beg, 'native'), y = y,
    width = unit(dgap$end - dgap$beg, 'native'), height = ht,
    just = c('left', 'center'),
    gp = gpar(lwd = 0.1, fill = fill, alpha = 0.9), vp = vp)
}
plot_link <- function(dcmp, dsnp, y = unit(0.5, 'npc'), height = 30, tgt.up = T,
  fill.p = 'skyblue1', fill.n = 'tomato', alpha = 0.5, vp = NULL) {
  comp = dcmp; snp = dsnp
  if(is.null(comp)) { return(NULL) }
  comp.xs = c()
  comp.fills = c()
  comp.ids = c()
  for (i in 1:nrow(comp) ) {
    if( comp$qsrd.a[i] == comp$tsrd.a[i] ) {
      comp.x = c(comp$tbeg.a[i], comp$tend.a[i], comp$qend.a[i], comp$qbeg.a[i])
      comp.fill = fill.p
      comp.id = rep(i, 4)
    } else {
      comp.x = c(comp$tbeg.a[i], comp$tend.a[i], comp$qbeg.a[i], comp$qend.a[i])
      comp.fill = fill.n
      comp.id = rep(i, 4)
    }
    comp.xs = c(comp.xs, comp.x)
    comp.fills = c(comp.fills, comp.fill)
    comp.ids = c(comp.ids, comp.id)
  }
  
  if(tgt.up) {
    link.y = c(1,1,-1,-1)
  } else {
    link.y = c(-1,-1,1,1)
  }
  tmp = rep(link.y, nrow(comp))
  
  linkgrob = polygonGrob(
    x = unit(comp.xs, 'native'),
    y = y + unit(tmp * height / 2, 'points'),
    id = comp.ids,
    gp = gpar(fill = comp.fills, alpha = alpha, lwd = 0, col = NA),
    vp = vp)
  snp = NULL ### do not plot SNPs
  if(is.null(snp)) {
    gList(linkgrob)
  } else {
  if(tgt.up) {
    tsnp.y0 = y + unit(rep(height/2, nrow(snp)), 'points')
    tsnp.y1 = y + unit(rep(height/2-1, nrow(snp)), 'points')
    qsnp.y0 = y - unit(rep(height/2, nrow(snp)), 'points')
    qsnp.y1 = y - unit(rep(height/2-1, nrow(snp)), 'points')
  } else {
    tsnp.y0 = y - unit(rep(height/2, nrow(snp)), 'points')
    tsnp.y1 = y - unit(rep(height/2-1, nrow(snp)), 'points')
    qsnp.y0 = y + unit(rep(height/2, nrow(snp)), 'points')
    qsnp.y1 = y + unit(rep(height/2-1, nrow(snp)), 'points')
  }
    gList(linkgrob,
  segmentsGrob(
    x0 = unit(snp$tpos.a, 'native'), x1 = unit(snp$tpos.a, 'native'), 
    y0 = tsnp.y0, y1 = tsnp.y1,
    gp = gpar(alpha = 1, lwd = 0.5), vp = vp),
  segmentsGrob(
    x0 = unit(snp$qpos.a, 'native'), x1 = unit(snp$qpos.a, 'native'), 
    y0 = qsnp.y0, y1 = qsnp.y1,
    gp = gpar(alpha = 1, lwd = 0.5), vp = vp)
    )
  }
}
plot_hist <- function(ds, col = 'grey', vp = NULL) {
  ds = ds[,c('beg.a', 'end.a', 'val')]
  colnames(ds)[1:3] = c("beg", "end", "score")
  segmentsGrob( 
    x0 = unit(ds$beg, 'native'), x1 = unit(ds$beg, 'native'), 
    y0 = 0, y1 = unit(ds$score, 'native'), 
    gp = gpar(lwd = 0.1, col = col, alpha = 1), vp = vp)
}
plot_hist_rect <- function(ds, col = 'grey20', vp = NULL) {
  ds = ds[,c('beg.a', 'end.a', 'val')]
  colnames(ds)[1:3] = c("beg", "end", "score")
  rectGrob( 
    x = unit(ds$beg, 'native'), width = unit(ds$beg-ds$end, 'native'), 
    y = 0, height = unit(ds$score, 'native'),  just = c('left', 'bottom'),
    gp = gpar(col = NA, fill = col, lwd = 0, alpha = 1), vp = vp)
}
plot_snp <- function(ds, y = unit(0.5, 'npc'), col = 'black', vp = NULL) {
  if(is.null(ds)) {return(NULL)}
  segmentsGrob( 
    x0 = unit(ds$pos.a, 'native'), x1 = unit(ds$pos.a, 'native'), 
    y0 = y - unit(1, 'points'), y1 = y + unit(1, 'points'), 
    gp = gpar(alpha = 1, lwd = 0.5), vp = vp)
}
plot_reads <- function(ds, col.p = 'navy', col.n = 'salmon', vp = NULL) {
  if(is.null(ds)) {return(NULL)}
  ds = ds[,c('beg.a', 'end.a', 'srd', 'y')]
  dsp = ds[ds$srd == '+',]
  dsn = ds[ds$srd == '-',]
  lwd = ifelse(max(ds$y) <= 15, 1, 0.8)
  
  grobs = gList()
  if(nrow(dsp) > 0) {
  grobs = gList(grobs,
  segmentsGrob( 
    x0 = unit(dsp$beg.a, 'native'), x1 = unit(dsp$end.a, 'native'), 
    y0 = unit(dsp$y, 'native'), y1 = unit(dsp$y, 'native'), 
    gp = gpar(lwd = lwd, col = col.p), vp = vp)
  )
  }
  if(nrow(dsn) > 0) {
  grobs = gList(grobs,
  segmentsGrob( 
    x0 = unit(dsn$beg.a, 'native'), x1 = unit(dsn$end.a, 'native'), 
    y0 = unit(dsn$y, 'native'), y1 = unit(dsn$y, 'native'), 
    gp = gpar(lwd = lwd, col = col.n), vp = vp)
  )
  }
  grobs
}
plot_scale <- function(max_len, x = unit(0.98, 'npc'), y = unit(0.2, 'npc'), largescale = F, vp = NULL) {
	if(largescale) {
  	len = diff( pretty(1:max_len, 40)[1:2] )
  	name = sprintf("%gMb", len / 1000000)
  } else {
  	len = diff( pretty(1:max_len, 20)[1:2] )
  	name = sprintf("%gKb", len / 1000)
  }
  scalegrob1 = segmentsGrob( 
    x0 = x - unit(len, 'native'), x1 = x, y0 = y, y1 = y, 
    vp = vp)
  scalegrob2 = segmentsGrob( 
    x0 = unit.c(x - unit(len, 'native'), x),
    x1 = unit.c(x - unit(len, 'native'), x),
    y0 = rep(y, 2), 
    y1 = rep(y + unit(3, 'points'), 2),
    vp = vp)
  scalegrob3 = textGrob( name,
    x = x - unit(len / 2, 'native'), 
    y = y + unit(5, 'points'), 
    just = c("center", "bottom"), 
    gp = gpar(cex = 0.75, fontfamily = "serif"), vp = vp)
  gList(scalegrob1, scalegrob2, scalegrob3)
}
plot_legend_gene <- function(x = unit(0.2, 'npc'), opt = 'all', vp = NULL) {
  stopifnot(opt %in% c("all", "nocrp", "nonbs"))
  if(opt == 'all') {
    fill = c('TE' = 'slategray3', 'Non-TE Gene' = 'brown4', 'NBS-LRR' = 'forestgreen', 'CRP' = 'dodgerblue')
  } else if(opt == 'nonbs') {
    fill = c('TE' = 'slategray3', 'Non-TE Gene' = 'brown4', 'CRP' = 'dodgerblue')
  } else {
    fill = c('TE' = 'slategray3', 'Non-TE Gene' = 'brown4', 'NBS-LRR' = 'forestgreen')
  }
  n = length(fill)
  ds = data.frame(fam = names(fill), col = as.character(fill), as.is = T)
  wds = nchar(as.character(ds$fam)) * 6
  xs = x + unit(c(1:n)*15 + c(0,cumsum(wds)[1:(n-1)]), 'points')
  rects = rectGrob( 
    x = xs, y = unit(0.5, 'npc'),
    width = unit(10, 'points'), height = unit(5, 'points'),
    just = c('left', 'center'),
    gp = gpar(lwd = 0, fill = as.character(ds$col), col = NA, alpha = 0.9), vp = vp)
  texts = textGrob( ds$fam, 
    x = xs + unit(15, 'points'), y = unit(0.5, 'npc'), 
    just = c("left", "center"), 
    gp = gpar(cex = 0.7, fontfamily = "serif"), vp = vp)
  gList(rects, texts)
}
plot_title <- function(main = '', subtitle = '', vp = NULL) {
  maingrob = textGrob(main, 
    x = unit(0.5, 'npc'), y = unit(0.7, 'npc'), 
    gp = gpar(fontface = 'bold', fontfamily = 'serif'),
    just = c("center", "top"), 
    vp = vp)
  subgrob = textGrob(subtitle, 
    x = unit(0.5, 'npc'), y = unit(0.7, 'npc') - unit(1, 'lines'), 
    gp = gpar(fontface = 'bold', fontfamily = 'serif'),
    just = c("center", "top"), 
    vp = vp)
  gList(maingrob, subgrob)
}


