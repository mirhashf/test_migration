// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.net.Socket;
import java.util.HashMap;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Responsible for sending the aligner output to the right node based on the genome location of the aligned reads.
 * TODO: Reads should be sent to sorters in larger packets.
 * @author Amirhossein Kiani (amir@binatechnologies.com) Created on Feb 13, 2012
 */
public class AlignmentDeMultiplexer extends Thread {
  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Thread#run()
   */

  private static final String NEW_LINE = System.getProperty("line.separator");
  private BufferedReader alignmentReader;
  private Map<GenomeBucket, Socket> genomeBucketMap;
  private static Logger logger = LoggerFactory.getLogger(AlignmentDeMultiplexer.class);
  private JobRequest jobRequest;

  public BufferedReader getReader() {
    return alignmentReader;
  }

  /**
   * 
   */
  public AlignmentDeMultiplexer() {
    genomeBucketMap = new HashMap<GenomeBucket, Socket>();
  }

  public void setBufferedReader(BufferedReader alignmentReader) {
    this.alignmentReader = alignmentReader;
  }


  public void addGenomeBucketMapping(GenomeBucket bucket, Socket socket) {
    genomeBucketMap.put(bucket, socket);
  }

  @Override
  public void run() {
    String line = null;
    try {
      while ((line = alignmentReader.readLine()) != null) {
        if (line.startsWith("@")) {
          logger.warn("Throwing away header line: " + line);
          continue;
        }
        for (GenomeBucket bucket : genomeBucketMap.keySet()) {
          if (bucket.contains(line)) {
            sendLineToSocket(line, genomeBucketMap.get(bucket));
          }
        }
      }

      alignmentReader.close();

      for (Socket socket : genomeBucketMap.values()) {
        socket.close();
      }

    } catch (IOException e) {
      logger.error("Failed to run alignment multiplexer.", e);
      logger.info("Stopping the job.");
      ZooKeeperService.killJobSilently(jobRequest.getId());
    }
  }


  private void sendLineToSocket(String read, Socket socket) throws IOException {
    logger.debug("Sending to " + socket + ": " + read);
    OutputStream os = socket.getOutputStream();
    os.write((read + NEW_LINE).getBytes());
    os.flush();
  }
}
