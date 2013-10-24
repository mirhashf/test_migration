function is_deletion(str1, str2,    length_diff, j) {
  length_diff = length(str2) - length(str1);
  if (length_diff <= 0) { return -1; }
  for (j=length(str1); j >= 0; j--) {
    if (substr(str1, 1, j) == substr(str2, 1, j) && substr(str1, j+1) == substr(str2, j+1+length_diff)) {
      return j;
    } 
  }
  return -1;
}

# For non-comment lines in the VCF
!/^#/{
  start = $2 - 1;
  filter = $7;

  ref_allele = toupper($4);
  alt_allele = toupper($5);
  event_type = "deletion"
  deletion_length = length(ref_allele) - length(alt_allele)
  deletion_offset = is_deletion(alt_allele, ref_allele)
  if (deletion_offset != -1) {
    start = start + deletion_offset
    deletion_count++;
  } else {
    non_deletion_count++
    non_deletion_records[non_deletion_count] = $0;
    event_type = "other"
    if (ignore_other == 1) {
      next;
    }
    #printf("### %s\t%d\n", $0, index(ref_allele, alt_allele));
    #next;
  }

  gt_field = $10;
  split(gt_field, gt_fields, ":");
  gt = gt_fields[1];
  is_phased = index(gt, "|");
  is_pass = index(gt_field, "PASS") > 0;
  if (is_phased) {
    phased_count++;
  } else {
    #next;
  }
  if (is_pass) {
    pass_count++;
  }
  if (filter == "0") {
    bad_filter++;
    bad_filter_records[bad_filter] = $0;
    #next;
  }
  end = start + deletion_length;
  svlen = deletion_length;
  split($8, info_split, ";");
  for (field in info_split) {
    if (match(info_split[field], "^END=")) {
      split(info_split[field], field_split, "=");
      end = field_split[2];
    }
    if (match(info_split[field], "^SVLEN=")) {
      split(info_split[field], field_split, "=");
      svlen = -field_split[2];
    }
  }
  if (svlen != deletion_length) {
    mismatching_svlen++;
  }
  #computed_svlen = start - end;
  if (end - start - deletion_length == 0) {
    end_match = "end_match"
  } else {
    end_match = "end_mismatch"
  }
  print $1 "\t" start "\t" end "\t" event_type "\t" 100.0 * (end - start - deletion_length) / deletion_length "\t" end_match "\t" deletion_offset "\t" gt "\t" $0 
  #print $1, start, end, svlen, "is_pass =", is_pass, computed_svlen, gt, $8, $10
  #if (computed_svlen != svlen) {
  #  mismatching_end++;
  #  mismatching_rows[mismatching_end] = $0;
  #}
}

END {
  #print "###==============================";
  #printf("### %d/%d svs deletions\n", deletion_count, NR);
  #printf("### %d/%d svs pass\n", pass_count, NR);
  #printf("### %d/%d svs phased\n", phased_count, NR);
  #printf("### %d/%d svs with filter value = 0\n", bad_filter, NR);
  #printf("### %d/%d svs with SVLEN not matching deletion_length\n", mismatching_svlen, NR);
  #printf("### %d/%d svs with end not matching start - svlen\n", mismatching_end, NR);
  #print "### Printing rows which are not deletions"
  #for (i in non_deletion_records) {
  #  print "###Nondeletion", non_deletion_records[i]
  #}
  #print "### Printing rows with filter column = 0"
  #for (i in bad_filter_records) {
  #  print "###", bad_filter_records[i];
  #}
  #print "### Printing rows with mismatching SVLEN, END, START"
  #for (i in mismatching_rows) {
  #  print "###", mismatching_rows[i];
  #}
}
