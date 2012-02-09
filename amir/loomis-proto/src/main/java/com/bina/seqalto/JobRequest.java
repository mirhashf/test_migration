// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

/**
 * @author Amirhossein Kiani (amir@binatechnologies.com) Created on Jan 29, 2012
 * 
 *         Represents a JSON {@link JobRequest} that comes from client for the leader. Currently for
 *         simplicity, every job should run on all the nodes.
 */
public class JobRequest {
  private String binary;
  private String input;
  private String output;
  private String id;
  private boolean started;


  /**
   * Constructor
   */
  public JobRequest() {
    setStarted(false);
    input = new String();
    output = new String();
    id = new String();
  }

  public boolean getStarted() {
    return started;
  }

  public void setStarted(boolean started) {
    this.started = started;
  }

  public String getBinary() {
    return binary;
  }

  public void setBinary(String binary) {
    this.binary = binary;
  }

  public String getInput() {
    return input;
  }

  public void setInput(String input) {
    this.input = input;
  }

  public String getOutput() {
    return output;
  }

  public void setOutput(String output) {
    this.output = output;
  }

  public String getId() {
    return id;
  }

  public void setId(String id) {
    this.id = id;
  }
}
