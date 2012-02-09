// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.io.File;
import java.io.IOException;

import org.apache.zookeeper.KeeperException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 * 
 *         Run's jobs on {@link JobRunner}s (in a bash shell for now...)
 */
public class JobRunner {
  private static Logger logger = LoggerFactory.getLogger(JobRunner.class);

  void run(JobRequest job, RunnerNode node) throws IOException, InterruptedException, KeeperException {
    String myDir = ZooKeeperService.getNodeName();

    logger.info("Running job " + job.getId() + " on node " + ZooKeeperService.getNodeName()
        + " in working dir " + new File(".").getCanonicalPath());
    // TODO: use netcat -l PORT > fifo to listen on a port


    Runtime.getRuntime().exec("mkdir -p " + myDir);
    String path = new File(".").getCanonicalPath() + "/" + myDir + "/";
    String[] cmd =
        {"/bin/sh", "-c",
            job.getBinary() + " 1>" + path + job.getId() + ".o" + " 2>" + path + job.getId() + ".e"};

    byte[] read = new byte[200];
    Runtime.getRuntime().exec(cmd).getInputStream().read(read);
    Thread.sleep(500);
    System.out.println(new String(read));
    ZooKeeperService.runnerDone(job);
  }

}
