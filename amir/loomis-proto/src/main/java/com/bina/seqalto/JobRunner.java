// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.Map;
import java.util.Map.Entry;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.util.StringUtils;

/**
 * Run's jobs. For now this supports only a WGS pipeline. Should make it more
 * general later to support other types of analyses.
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 */
public class JobRunner extends Thread {
  private static Logger logger = LoggerFactory.getLogger(JobRunner.class);
  private JobRequest job;
  private RunnerNode node;
  private ServerSocket leaderSocket, sorterSocket;
  private static final String NEW_LINE = System.getProperty("line.separator");
  private Sorter sorter;

  public JobRunner(JobRequest job, RunnerNode node) {
    this.job = job;
    this.node = node;
  }

  @Override
  public void run() {
    BufferedReader inputData;
    String inputLine;

    try {
      if (job.isStreaming()) {
        leaderSocket = NetworkUtils.openServerSocketWithRetry(node.getPort());
        sorterSocket = NetworkUtils.openServerSocketWithRetry(node.getSorterPort());

        // create and initialize sorter
        sorter = new Sorter(job);
        sorter.setServerSocket(sorterSocket);
        sorter
            .setGenomeBucket(getMyOwnGenomeBucket(node, node.getGenomeBucketToRunnerNodeMapping()));
        sorter.setNumberOfAligners(node.getGenomeBucketToRunnerNodeMapping().keySet().size());
        sorter.start();

        logger.debug("Successfully bound to socket on port " + node.getPort());
        leaderSocket.setReuseAddress(true);
        sorterSocket.setReuseAddress(true);

        node.setListening(true);
        ZooKeeperService.setRunnerNode(node, job.getId());

        Socket clientSocket = leaderSocket.accept();
        inputData = new BufferedReader(new InputStreamReader(clientSocket.getInputStream()));
      } else {
        inputData =
            new BufferedReader(
                new FileReader(new File(job.getProperty(JobParameters.FASTQ_PAIR_1))));
      }

      logger.info("Running job " + job.getId() + " on node " + ZooKeeperService.getNodeName()
          + " in working dir " + new File(".").getCanonicalPath());

      Runtime.getRuntime().exec("mkdir -p " + job.getProperty(JobParameters.OUTPUT_DIR));

      String[] cmd = {"/bin/sh", "-c", SpringPropertiesUtil.getProperty("seqalto.binary")
      // + " -idx " + job.getProperty(JobParameters.SEQALTO_INDEX)

          };

      logger.info("Running " + StringUtils.arrayToDelimitedString(cmd, " "));

      Process process = Runtime.getRuntime().exec(cmd);
      OutputStream processInput = process.getOutputStream();


      logger.info("Starting an aligner multiplexer to multiplex the output.");
      // mutliplexer to send the right read to the right place
      AlignmentDeMultiplexer multiplexer = new AlignmentDeMultiplexer();
      multiplexer.setBufferedReader(new BufferedReader(new InputStreamReader(process
          .getInputStream())));

      for (Entry<GenomeBucket, String> entry : node.getGenomeBucketToRunnerNodeMapping().entrySet()) {
        GenomeBucket bucket = entry.getKey();
        String runnerNodeId = entry.getValue();
        RunnerNode runnerNode = ZooKeeperService.getRunnerNode(runnerNodeId, job.getId());
        Socket nodeSocket =
            NetworkUtils.openSocketWithRetry(runnerNode.getHost(), runnerNode.getSorterPort());
        multiplexer.addGenomeBucketMapping(bucket, nodeSocket);
      }

      // start listening on multiplexer
      multiplexer.start();



      logger.debug("Starting to read data from input source...");

      while ((inputLine = inputData.readLine()) != null && !Thread.interrupted()) {
        logger.debug("From Leader: " + inputLine);
        processInput.write((inputLine + NEW_LINE).getBytes());
        processInput.flush();
      }

      if (Thread.interrupted()) {
        logger.info("Job has been killed! Intrrupting the sorter and destroying the process.");
        sorter.interrupt();
        process.destroy();
      }

      processInput.close();

      // wait for alignment to be done
      process.waitFor();


      logger.debug("Waiting for alignment multiplexer to be done.");
      sorter.join();

      logger.debug("Waiting for sorter to be done.");
      multiplexer.join();

      logger.debug("Closing input socket.");
      leaderSocket.close();

      logger.debug("Closing sorter's socket.");
      sorterSocket.close();
      ZooKeeperService.runnerDone(job);
    } catch (Exception e) {
      logger.error("Job " + job.getId() + " faild!", e);
      // tell everyone to die!
      ZooKeeperService.killJobSilently(job.getId());
    }
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#finalize()
   */
  @Override
  protected void finalize() throws Throwable {
    super.finalize();
    leaderSocket.close();
  }


  private GenomeBucket getMyOwnGenomeBucket(RunnerNode node, Map<GenomeBucket, String> mapping) {
    for (Entry<GenomeBucket, String> entry : mapping.entrySet()) {
      if (entry.getValue().equals(node.getId())) {
        return entry.getKey();
      }
    }
    return null;
  }
}
