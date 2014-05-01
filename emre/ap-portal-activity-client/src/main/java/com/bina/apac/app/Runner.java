package com.bina.apac.app;

import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.webapp.WebAppContext;

/**
 * This class is provided for convenience when 
 * running the app from an IDE like Eclipse
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
