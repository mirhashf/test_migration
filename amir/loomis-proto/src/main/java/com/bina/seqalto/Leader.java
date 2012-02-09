// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import javax.inject.Inject;

import org.springframework.stereotype.Service;

/**
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 *
 * Leader's work is done in this class! Clients send {@link JobRequest}s to Leader if they receive one.
 */
@Service
public class Leader {
  @Inject
  ZooKeeperService zooKeeper;
  public synchronized void process(JobRequest request) throws Exception{
    // create job node
    zooKeeper.createJob(request);
  }
}
