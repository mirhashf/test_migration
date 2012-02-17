// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.Socket;
import java.util.List;

import javax.inject.Inject;

import org.apache.zookeeper.KeeperException;
import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Service;

/**
 * Leader's work is done in this class. Currently responsible for: \n1- Streaming input to followers
 * for streaming jobs by spawning a {@link InputStreamer} thread. \n2- Killing jobs through (only
 * leader can kill a job if it's accessed externally)
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 */
@Service
public class Leader {
  @Inject
  ZooKeeperService zooKeeper;

  private static Logger logger = LoggerFactory.getLogger(Leader.class);
  private static final String NEW_LINE = System.getProperty("line.separator");
  private InputStreamer inputStreamer;

  public InputStreamer getInputStreamer() {
    return inputStreamer;
  }

  public void setInputStreamer(InputStreamer inputStreamer) {
    this.inputStreamer = inputStreamer;
  }

  public void process(JobRequest request) throws Exception {
    // create job node
    zooKeeper.createJob(request);
  }

  public void streamInput(JobRequest request) throws Exception {
    inputStreamer = new InputStreamer();
    inputStreamer.setRequest(request);
    logger.debug("Starting streamer thread from leader for job " + request.getId());
    inputStreamer.start();
  }


  /**
   * Thread to stream input to followers from leader.
   * 
   * @author Amirhossein Kiani (amir@binatechnologies.com) Created on Feb 16, 2012
   */
  class InputStreamer extends Thread {
    private JobRequest request;

    public JobRequest getRequest() {
      return request;
    }

    public void setRequest(JobRequest request) {
      this.request = request;
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.lang.Thread#run()
     */
    @Override
    public void run() {
      try {
        logger.info("Leader started input streaming thread for job " + request.getId());
        int packetSizeInLines = request.getPacketSizeInLines();
        List<RunnerNode> runners = ZooKeeperService.getJobRunners(request);
        File inputFile = new File(request.getProperty(JobParameters.FASTQ_PAIR_1));
        BufferedReader reader = new BufferedReader(new FileReader(inputFile));
        StringBuilder packet = new StringBuilder();
        String line;
        int lineNumber = 1;
        long packetNumber = 0;

        int numRunners = runners.size();

        Socket[] runnerSockets = new Socket[numRunners];

        int i = 0;
        for (RunnerNode runner : runners) {
          runnerSockets[i] = NetworkUtils.openSocketWithRetry(runner.getHost(), runner.getPort());
          i++;

        }

        PrintWriter out = null;

        while ((line = reader.readLine()) != null && !Thread.interrupted()) {
          lineNumber++;
          packet.append(line + NEW_LINE);
          // stream packets to runners one by one.
          if (lineNumber % packetSizeInLines == 0) {
            int runnerIndex = (int) (packetNumber % numRunners);
            RunnerNode runner = runners.get(runnerIndex);

            logger.debug("Sending packet " + packetNumber + " to runner " + runner.getHost()
                + " port " + runner.getPort() + "...");

            // actually send the packet here:
            out = new PrintWriter(runnerSockets[runnerIndex].getOutputStream(), true);
            out.write(packet.toString());
            out.flush();

            // clear the "buffer"
            packetNumber++;
            packet = new StringBuilder();
          }
        }

        // close all the sockets when we are done
        for (i = 0; i < numRunners; i++) {
          runnerSockets[i].close();
          logger.debug("Closed the socket to runner " + runners.get(i).getHost() + " on port "
              + runners.get(i).getPort());
        }
      } catch (Exception e) {

        logger.error("Error in streaming input to runners!.. killing job " + request.getId(), e);
        ZooKeeperService.killJobSilently(request.getId());

      }
    }
  }


  /**
   * Kill a job
   * 
   * @param jobId
   * @throws IOException
   * @throws InterruptedException
   * @throws KeeperException
   * @throws JsonMappingException
   * @throws JsonParseException
   */
  public void killJob(String jobId) throws JsonParseException, JsonMappingException,
      KeeperException, InterruptedException, IOException {
    // oh man... this is going to be challenging...
    // this will tell all the runners to stop doing whatever they are doing..
    ZooKeeperService.killJob(jobId);
  }


}
