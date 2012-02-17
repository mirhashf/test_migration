// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.util.HashMap;
import java.util.Map;

/**
 * Data structure used to keep the configuration for a runner (Follower) in the ZooKeeper tree.
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com) Created on Jan 29, 2012
 */
public class RunnerNode {
  private int port;
  private int sorterPort;

  public int getSorterPort() {
    return sorterPort;
  }

  public void setSorterPort(int sorterPort) {
    this.sorterPort = sorterPort;
  }

  private String host;
  private boolean started = false;
  private String id;
  private Map<String, String> properties;
  private boolean listening;


  private Map<GenomeBucket, String> genomeBucketToRunnerNodeMapping =
      new HashMap<GenomeBucket, String>();

  public boolean isListening() {
    return listening;
  }

  public void setListening(boolean listening) {
    this.listening = listening;
  }

  public RunnerNode() {
    properties = new HashMap<String, String>();
  }

  public void addProperty(String key, String value) {
    properties.put(key, value);
  }

  public void getProperty(String key) {
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

  public int getPort() {
    return port;
  }

  public void setPort(int l) {
    this.port = l;
  }

  public Map<GenomeBucket, String> getGenomeBucketToRunnerNodeMapping() {
    return genomeBucketToRunnerNodeMapping;
  }

  public void setGenomeBucketToRunnerNodeMapping(
      Map<GenomeBucket, String> genomeBucketToRunnerNodeMapping) {
    this.genomeBucketToRunnerNodeMapping = genomeBucketToRunnerNodeMapping;
  }

  public void addGenomeBucketToRunnerNodeMapping(GenomeBucket bucket, String runnerNode) {
    this.genomeBucketToRunnerNodeMapping.put(bucket, runnerNode);
  }

  public void removeGenomeBucket(GenomeBucket bucket) {
    this.genomeBucketToRunnerNodeMapping.remove(bucket);
  }



}
