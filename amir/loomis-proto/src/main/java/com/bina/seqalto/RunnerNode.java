// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

/**
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 * Created on Jan 29, 2012
 * 
 * Represents a node on the box
 */
public class RunnerNode {
  private long ip;
  private String host;
  
  public String getHost() {
    return host;
  }
  public void setHost(String host) {
    this.host = host;
  }
  public long getIp() {
    return ip;
  }
  public void setIp(long l) {
    this.ip = l;
  }
}
