require(Biostrings)

aa_pw_dist <- function(rw) {
  qseq = rw['qseq']
  tseq = rw['tseq']
  qlen = as.numeric(nchar(qseq))
  tlen = as.numeric(nchar(tseq))
  pw = pairwiseAlignment(AAString(tseq), AAString(qseq), type = "global",
    substitutionMatrix = "BLOSUM62", gapOpening = -3, gapExtension = -1)
  mat = nmatch(pw)
  mis = nmismatch(pw)
  indel = nindel(pw)
  qgap = indel@insertion[2]
  tgap =indel@deletion[2]
  alen = nchar(pw)
  stopifnot(alen == mis + mat + qgap + tgap)
  qres = qlen - (mat + mis + tgap)
  tres = tlen - (mat + mis + qgap)
  qgap = qgap + tres
  tgap = tgap + qres
  len = alen + qres + tres
  stopifnot(len == mis + mat + qgap + tgap)
  c('alen'=alen,'mat'=mat,'mis'=mis,'tgap'=tgap,'qgap'=qgap,'len'=len)
}
