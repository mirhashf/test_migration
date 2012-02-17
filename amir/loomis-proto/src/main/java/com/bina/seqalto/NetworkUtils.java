// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.net.ConnectException;
import java.net.ServerSocket;
import java.net.Socket;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Utility class for common network related methods
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com) Created on Feb 13, 2012
 */
public class NetworkUtils {
  private static Logger logger = LoggerFactory.getLogger(NetworkUtils.class);

  public final static int NUMBER_OF_CONNECTION_RETIES = 5;
  public final static int WAIT_TIME_BETWEEN_CONNECTION_RETIES = 10000;

  public static Socket openSocketWithRetry(String host, int port) throws Exception {
    Socket socket = null;
    for (int i = 1; i <= NUMBER_OF_CONNECTION_RETIES; i++) {
      try {
        logger.debug("Openning client socket to write output... (attempt " + i + ")");
        socket = new Socket(host, port);
        return socket;
      } catch (ConnectException ex) {
        logger.debug("Connection refused. Trying again in " + WAIT_TIME_BETWEEN_CONNECTION_RETIES
            + " miliseconds...");
        Thread.sleep(WAIT_TIME_BETWEEN_CONNECTION_RETIES);
      }
    }
    throw new ConnectException("Ran out of" + "(" + NUMBER_OF_CONNECTION_RETIES + ") retries!...");
  }

  public static ServerSocket openServerSocketWithRetry(int port) throws Exception {
    ServerSocket serverSocket;
    for (int i = 1; i <= NetworkUtils.NUMBER_OF_CONNECTION_RETIES; i++) {
      try {
        logger.debug("Openning server socket to read input... (attempt " + i + ")");
        serverSocket = new ServerSocket(port);
        return serverSocket;
      } catch (ConnectException ex) {
        logger.debug("Connection refused. Trying again in "
            + NetworkUtils.WAIT_TIME_BETWEEN_CONNECTION_RETIES + " miliseconds...");
        Thread.sleep(NetworkUtils.WAIT_TIME_BETWEEN_CONNECTION_RETIES);
      }
    }
    throw new ConnectException("Ran out of" + "(" + NUMBER_OF_CONNECTION_RETIES + ") retries!...");
  }
}
