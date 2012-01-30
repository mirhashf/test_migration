// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;


/**
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 *
 * Leader's work is done in this class! Clients send {@link JobRequest}s to Leader if they receive one.
 */
public class Leader {
  public void process(JobRequest request) throws Exception{
    // create job node
    ZooKeeperService.createJob(request);

  }
}
