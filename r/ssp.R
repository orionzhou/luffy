read_ssp <- function(fi, quiet=T) {
  line1 = scan(fi, what=integer(0), sep="\t", nlines=1, quiet=quiet)
  n_ind = line1[1]
  n_pos = line1[2]
  poss = scan(fi, what=double(0), sep="\t", nlines=1, skip=1, quiet=quiet)
  if(n_pos != length(poss)) {
    warning(paste("line 2 not", n_pos, "positions:", length(poss)))
  }
  anc = scan(fi, what=character(0), sep="\t", nlines=1, skip=2, quiet=quiet)
  if(n_pos != length(anc)) {
    warning(paste("line 3 not", n_pos, "ancestral alleles:", length(anc)))
  }
  hasOG = sum( anc %in% c('N', '?') ) < n_pos
  inds = c()
  data = matrix(NA, nrow=n_ind, ncol=n_pos)
  for (i in 1:n_ind) {
    tmp = scan(fi, what=character(), sep="\t", nlines=1, skip=i+2, quiet=quiet)
    inds = c(inds, tmp[1])
    data[i,] = tmp[-1]
  }
  rownames(data) = inds
  if(hasOG) {
    n_ind = n_ind + 1
    inds = c('anc', inds)
    data = rbind(anc, data)
  }
  list(n_ind=n_ind, n_pos=n_pos, hasOG=hasOG, inds=inds, poss=poss, data=data)
}
get_sfs_ssp = function(ssp) {
  sfs = matrix(NA, ncol=4, nrow=ssp$n_pos)
  colnames(sfs) = c("n_N", "n_anc", "n_der", "freq_der")
  rownames(sfs) = ssp$poss
  for (i in 1:ssp$n_pos) {
    if(ssp$hasOG) {
      anc = ssp$data[1,i]
      alleles = unique(ssp$data[,i])
      t = table(ssp$data[2:ssp$n_ind,i])
      if(anc == 'N') {
        warning(paste("ancestral state of position", ssp$poss[i], "is N"))
      } else if(sum(! alleles %in% c(anc, 'N')) != 1) {
        warning(paste("position", ssp$poss[i], "not having 2 alleles"))
      } else {
        der = alleles[which(!alleles %in% c(anc,'N'))]
        if(sum(alleles=='N')>0) sfs[i,'n_N'] = t['N'] else sfs[i, 'n_N'] = 0
        sfs[i, 'n_anc'] = t[anc]
        sfs[i, 'n_der'] = t[der]
      }
    } else {
      t = table(ssp$data[,i])
      alleles = names(sort(t)) 
      if(sum(alleles != 'N') != 2) {
        warning(paste("position", ssp$poss[i], "not having 2 alleles"))
      } else {
        alleles_not_N = alleles[which(alleles != 'N')]
        anc = alleles_not_N[2]
        der = alleles_not_N[1]
        if(sum(alleles=='N')>0) sfs[i,'n_N'] = t['N'] else sfs[i,'n_N'] = 0
        sfs[i, 'n_anc'] = t[anc]
        sfs[i, 'n_der'] = t[der]
      }
    }
  }
  sfs[,'freq_der'] = sfs[,'n_der'] / (sfs[,'n_anc'] + sfs[,'n_der'])
  sfs
}
get_one_sfs = function(nts, hasOG) {
  i = 0
  if(hasOG) {
    anc = nts[1]
    alleles = as.numeric(nts[-1])
    n_ind = length(nts) - 1
    if(anc == 0) {
      stop(paste("ancestral state of snp", i, "is N"))
    }
  } else {
    alleles = as.numeric(nts)
    n_ind = length(nts)
    t = sort(table(alleles[alleles!=0]))
    anc = names(t)[-1]
  }
  n_N = sum(alleles == 0)
  states = unique(alleles[alleles!=0])
  if(length(states) != 2) {
    stop(paste("snp", i, "does not have 2 alleles"))
  }
  der = states[states!=anc]
  n_anc = sum(alleles == anc)
  n_der = sum(alleles == der) 
  c(n_N, n_anc, n_der)
}

testHWE <- function(AA,Aa,aa) {
 total = AA+Aa+aa
 PA = (AA+Aa/2)/total
 Pa = 1-PA
 AAe = total * PA^2
 Aae = total * PA * Pa * 2
 aae = total * Pa * Pa
 chsq = (AA-AAe)^2/AAe + (Aa-Aae)^2/Aae + (aa-aae)^2/aae
 1-pchisq(chsq,df=1)
}
