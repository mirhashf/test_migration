// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.io.IOException;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.annotate.JsonCreator;
import org.codehaus.jackson.annotate.JsonProperty;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;

/**
 * Represents a genomic region (reference, start, end)
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com) Created on Feb 8, 2012
 */
public class GenomeRegion implements Comparable<GenomeRegion> {
  @JsonProperty
  private String reference;
  @JsonProperty
  private long start;
  @JsonProperty
  private long end;


  /**
   * Cosnstructor
   */
  public GenomeRegion() {
    // Had to have this to be able to conver the objet to JSON...
  }

  /**
   * Create a new Genomic Region
   * 
   * @param reference
   * @param start
   * @param end
   */
  public GenomeRegion(String reference, long start, long end) {
    super();
    this.reference = reference;
    this.start = start;
    this.end = end;
  }

  public String getReference() {
    return reference;
  }

  public void setReference(String reference) {
    this.reference = reference;
  }

  public long getStart() {
    return start;
  }

  public void setStart(long start) {
    this.start = start;
  }

  public long getEnd() {
    return end;
  }

  public void setEnd(long end) {
    this.end = end;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Comparable#compareTo(java.lang.Object)
   */
  public int compareTo(GenomeRegion o) {
    return Long.valueOf(this.end - this.start).compareTo(o.end - o.start);
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  public String toString() {
    return "{\"reference\":\"" + this.getReference() + "\", \"start\":" + this.getStart()
        + ", \"end\":" + this.getEnd() + "}";
  }

  public long getSize() {
    return this.end - this.start + 1;
  }

  /**
   * Analyzes a SAM record in ASCI format and returns true if this region contains the input read.
   * 
   * @param recordSamString
   * @return
   * @throws IOException
   */
  public boolean contains(String recordSamString) throws IOException {
    if (recordSamString.startsWith("@")) {
      // header... shouldn't have received this
      throw new IOException("SAM record text is invalid: " + recordSamString);
    }
    try {
      String[] record = recordSamString.split("\t", 5);
      String reference = record[2];
      long startLocation = Long.valueOf(record[3]);
      if (reference.equals(this.reference) && startLocation >= getStart()
          && startLocation < getEnd()) {
        return true;
      } else {
        return false;
      }
    } catch (Exception ex) {
      throw new IOException("SAM record text is invalid: " + recordSamString, ex);
    }
  }

  @JsonCreator
  public static GenomeRegion fromJSON(String val) throws JsonParseException, JsonMappingException,
      IOException {
    ObjectMapper mapper = new ObjectMapper();
    GenomeRegion region = mapper.readValue(val, GenomeRegion.class);
    return region;
  }

}
