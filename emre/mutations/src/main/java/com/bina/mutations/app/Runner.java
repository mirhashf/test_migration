package com.bina.mutations.app;

import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.webapp.WebAppContext;

/**
 * Helps with starting from within an IDE 
 * 
 * */
public class Runner {

  public static void main(String[] args) throws Exception {
    final Server server = new Server(9000);
    final WebAppContext context = new WebAppContext();
    context.setDescriptor("src/main/resources/webapp/WEB-INF/web.xml");
    context.setResourceBase("src/main/webapp");
    context.setContextPath("/");
    context.setParentLoaderPriority(true);

    server.setHandler(context);
    server.start();
    server.join();
  }
}
