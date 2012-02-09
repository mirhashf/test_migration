// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 * Created on Jan 29, 2012
 * 
 * Represents a node on the box
 */
public class RunnerNode {
  private long ip;
  private String host;
  private boolean started = false;
  private String id;
  private Map<String, String> properties;
  
  public RunnerNode() {
    properties  = new HashMap<String, String>();
  }
  
  public void addProperty(String key, String value){
    properties.put(key, value);
  }
  
  public void getProperty(String key){
    properties.get(key);
  }
  
  public String getId() {
    return id;
  }
  public void setId(String id) {
    this.id = id;
  }
  public boolean isStarted() {
    return started;
  }
  public void setStarted(boolean started) {
    this.started = started;
  }
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
