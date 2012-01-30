// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

/**
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 * Created on Jan 29, 2012
 * 
 * Represents a JSON {@link JobRequest} that comes from client for the leader.
 */
public class JobRequest {
  private String binary;
  private String input;
  private String output;
  private int id;
 
  
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

  public int getId() {
    return id;
  }

  public void setId(int id) {
    this.id = id;
  }
}
