// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.codehaus.jackson.annotate.JsonProperty;

/**
 * Represents a group of {@link GenomeRegion}'s.
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com) Created on Feb 9, 2012
 */

public class GenomeBucket {
  @JsonProperty
  List<GenomeRegion> regions;
  @JsonProperty
  long totalSize;

  public long getTotalSize() {
    return totalSize;
  }

  public void setTotalSize(long totalSize) {
    this.totalSize = totalSize;
  }

  /**
   * @param regions
   * @param totalSize
   */
  public GenomeBucket(List<GenomeRegion> regions, long totalSize) {
    super();
    this.regions = regions;
    this.totalSize = totalSize;
  }

  /**
   * Constructor
   */
  public GenomeBucket() {
    regions = new ArrayList<GenomeRegion>();
  }

  public void addRegion(GenomeRegion region) {
    totalSize += region.getSize();
    regions.add(region);
  }

  public List<GenomeRegion> getRegions() {
    return regions;
  }

  public long getSize() {
    return totalSize;
  }

  /**
   * Analyzes a SAM record in ASCI format and returns true if this region contains the input read.
   * 
   * @param samRecordString
   * @return
   * @throws IOException
   */
  public boolean contains(String samRecordString) throws IOException {
    for (GenomeRegion region : getRegions()) {
      if (region.contains(samRecordString)) {
        return true;
      }
    }
    return false;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  public String toString() {
    return "{\"regions\":" + regions.toString() + ",\"totalSize\":" + totalSize + "}";
  }

  
}
