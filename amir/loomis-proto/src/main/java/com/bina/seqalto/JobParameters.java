// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.util.HashMap;
import java.util.Map;

/**
 * Enum for parameters used in a {@link JobRequest}'s properties map.
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com) Created on Feb 13, 2012
 */
public enum JobParameters {

  FASTA_FILE("-fasta"), OUTPUT_DIR("-output-dir"), FASTQ_PAIR_1("-fastq_1"), FASTQ_PAIR_2("-fastq_2"), SEQALTO_INDEX("-seqalto-index");

  private final String flag;
  private static final Map<String, JobParameters> lookup = new HashMap<String, JobParameters>();
  static {
    for (JobParameters d : JobParameters.values())
      lookup.put(d.getFlag(), d);
  }

  private JobParameters(String flag) {
    this.flag = flag;
  }

  public String getFlag() {
    return flag;
  }

  public static JobParameters get(String flag) {
    return lookup.get(flag);
  }

}
