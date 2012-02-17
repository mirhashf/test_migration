// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Amirhossein Kiani (amir@binatechnologies.com) Created on Jan 29, 2012
 * 
 *         Represents a job that needs to run. It is also used to store a job's data structure on ZooKeeper.
 */
public class JobRequest {
  private String id;
  private boolean started;
  private boolean streaming;
  private boolean streamingStarted;
  private boolean stopped;
  
  public boolean isStreamingStarted() {
    return streamingStarted;
  }

  public void setStreamingStarted(boolean streamingStarted) {
    this.streamingStarted = streamingStarted;
  }

  private int packetSizeInLines = 1;
  private Map<JobParameters, String> jobProperties = new HashMap<JobParameters,String>();

  public Map<JobParameters, String> getJobProperties() {
    return jobProperties;
  }

  public void setJobProperties(Map<JobParameters, String> jobProperties) {
    this.jobProperties = jobProperties;
  }

  public int getPacketSizeInLines() {
    return packetSizeInLines;
  }

  public void setPacketSizeInLines(int packetSizeInLines) {
    this.packetSizeInLines = packetSizeInLines;
  }

  public boolean isStreaming() {
    return streaming;
  }

  public void setStreaming(boolean streaming) {
    this.streaming = streaming;
  }

  /**
   * Constructor
   */
  public JobRequest() {
    setStarted(false);
    id = new String();
  }

  public boolean getStarted() {
    return started;
  }

  public void setStarted(boolean started) {
    this.started = started;
  }
  public String getId() {
    return id;
  }

  public void setId(String id) {
    this.id = id;
  }
  
  public void setProperty(JobParameters property, String value){
    this.jobProperties.put(property, value);
  }
  
  public String getProperty(JobParameters property){
    return this.jobProperties.get(property);
  }

  public boolean isStopped() {
    return stopped;
  }

  public void setStopped(boolean stopped) {
    this.stopped = stopped;
  }
}
