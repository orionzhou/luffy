#!/bin/awk -f
function join( array, start, end, sep, result, i) {
  if (seq == "")
    sep = " ";
  else if (sep == SUBSEP)
    sep = "";
  result = array[start];
  for (i = start + 1; i <= end; i++) 
    result = result sep array[i];
  return result;
}
BEGIN {
  nPos = 10000;
}
{
  if(NR == 1) {
    print $1 " " nPos;
  } else if (NR == 2 || NR == 3) {
    for (i = 1; i<= nPos; i++) {
      a[i] = $i;
    }
    print join(a, 1, nPos, " ");
  } else {
    for (i = 1; i<= nPos + 1; i++) {
      a[i] = $i;
    }
    print join(a, 1, nPos+1, " ", result);
  }
}
