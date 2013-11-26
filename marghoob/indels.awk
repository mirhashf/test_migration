function is_deletion(str1, str2) {
  length_diff = length(str2) - length(str1);
  if (length_diff <= 0) { return -1; }
  for (j=0; j <= length(str1); j++) {
    if (substr(str1, 1, j) == substr(str2, 1, j) && substr(str1, j+1) == substr(str2, j+1+length_diff)) {
      return j;
    } 
  }
  return -1;
}

BEGIN {
  if (length(min_size) == 0) { min_size = 0 }
}

!/^#/ {
  ref_allele = $4
  alt_alleles = $5
  event = ""
  gt = $10
  split(gt, gt_split, ":")

  if (index(alt_alleles, ",")) {
    split(alt_alleles, alt_allele, ",");
  } else {
    alt_allele[1] = alt_alleles;
  }

  for (i in alt_allele) {
    begin = $2; end = $2 + 1;
    if (is_deletion(alt_allele[i], ref_allele) != -1) {
      events[i] = "deletion";
      begin = $2 + is_deletion(alt_allele[i], ref_allele) - 1;
      end = begin + length(ref_allele) - length(alt_allele[i]);
    } else if (is_deletion(ref_allele, alt_allele[i]) != -1) {
      begin = $2 + is_deletion(ref_allele, alt_allele[i]) - 1;
      end = begin + 1;
      events[i] = "insertion";
    } else if (length(alt_allele[i]) == 1 && length(ref_allele) == 1) {
      events[i] = "SNP";
    } else {
      end = begin + length(ref_allele) - length(alt_allele[i])
      events[i] = "other";
    }
    if (events[i] != "SNP" && (end - begin >= min_size)) {
      print $1"\t"begin"\t"end"\t"events[i]
    }
    #print $1"\t"begin"\t"end"\t"NR"\t"i","events[i]","(length(alt_allele[i]) - length(ref_allele))"\t"gt_split[1] "\t" $0
  }
  delete alt_allele
}
