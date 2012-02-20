// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.ArrayList;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.util.StringUtils;


/**
 * Wrapper thread responsible for running the sorter binary. Acts as a server and spawns a new
 * {@link SorterListener} thread for each aligner's request (up to the number of runner nodes) It
 * also waits for all the runner nodes to finish so it acts as a barrier for the sorters.
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com) Created on Feb 13, 2012
 */
public class Sorter extends Thread {

  private ServerSocket serverSocket;
  private static final Logger logger = LoggerFactory.getLogger(Sorter.class);
  private GenomeBucket genomeBucket = new GenomeBucket();
  private int numberOfConnectedAligners = 0;
  public static final String END_ALIGNMENT_SIGNAL = "END_ALIGNMENT";
  private int numberOfAligners = 1;
  private OutputStream sorterOutputStream;
  private static final String NEW_LINE = System.getProperty("line.separator");
  private JobRequest jobRequest;
  private ArrayList<SorterListener> listeners = new ArrayList<SorterListener>();

  /**
   * 
   */
  public Sorter(JobRequest jobRequest) {
    this.jobRequest = jobRequest;
  }

  public int getNumberOfAligners() {
    return numberOfAligners;
  }

  public void setNumberOfAligners(int numberOfAligners) {
    this.numberOfAligners = numberOfAligners;
  }

  public GenomeBucket getGenomeBucket() {
    return genomeBucket;
  }

  public void setGenomeBucket(GenomeBucket genomeBucket) {
    this.genomeBucket = genomeBucket;
  }

  public ServerSocket getServerSocket() {
    return serverSocket;
  }

  public void setServerSocket(ServerSocket inputSocket) {
    this.serverSocket = inputSocket;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Thread#run()
   */
  @Override
  public void run() {
    try {
      logger.info("Starting sorter on " + serverSocket);
      String sorterBinary = SpringPropertiesUtil.getProperty("sorter.binary");
      String[] cmd = {"/bin/sh", "-c", sorterBinary};
      Process process;
      process = Runtime.getRuntime().exec(cmd);

      logger.info("Running " + StringUtils.arrayToDelimitedString(cmd, " "));

      sorterOutputStream = process.getOutputStream();
      logger.info("Sorting on Genome Bucket " + genomeBucket);

      while (numberOfConnectedAligners < numberOfAligners) {
        numberOfConnectedAligners++;
        Socket socket = serverSocket.accept();
        logger.info("Starting Sorter Listener socket. " + numberOfConnectedAligners
            + " aligners out of " + numberOfAligners + "  have connected so far.");
        SorterListener listener = new SorterListener(jobRequest);
        listener.setSocket(socket);
        listener.start();
        listeners.add(listener);
      }

      logger.info("All aligners have been connected to the sorters.");

      logger.debug("Waiting for all sorter listener threads to be done...");
      for (SorterListener listener : listeners) {
        listener.join();
      }
      logger.debug("All sorter listener threads are done!");



    } catch (Exception e) {
      logger.error(
          "Failed to create sorter on socket (this can be ignored if job has been killed.)"
              + serverSocket, e);
      ZooKeeperService.killJobSilently(jobRequest.getId());
    }
  }

  private void writeToSorterProcess(String line) throws IOException {
    sorterOutputStream.write(line.getBytes());
    sorterOutputStream.flush();
  }


  /**
   * Listener thread for the sorter server that listens to the connected aligners.
   * 
   * @author Amirhossein Kiani (amir@binatechnologies.com) Created on Feb 16, 2012
   */
  private class SorterListener extends Thread {
    Socket socket;
    JobRequest jobRequest;

    public SorterListener(JobRequest jobRequest) {
      this.jobRequest = jobRequest;
    }

    public void setSocket(Socket socket) {
      this.socket = socket;
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.lang.Thread#run()
     */
    @Override
    public void run() {
      try {
        BufferedReader reader = new BufferedReader(new InputStreamReader(socket.getInputStream()));
        String line;
        while ((line = reader.readLine()) != null) {
          logger.debug("Recieved this in sorter listener: " + line);
          writeToSorterProcess(line + NEW_LINE);
        }

        logger.debug("Closing sorter listener socket.");
        socket.close();
      } catch (Exception e) {
        logger.error("Failed to create sorter listener on socket " + socket, e);
        ZooKeeperService.killJobSilently(jobRequest.getId());
      }
    }
  }
}
