package com.bina.apac.config;

import org.glassfish.jersey.jackson.JacksonFeature;
import org.glassfish.jersey.server.ResourceConfig;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ApplicationConfig extends ResourceConfig {

  private static final Logger logger = LoggerFactory.getLogger(ApplicationConfig.class);

  public ApplicationConfig() {
    logger.info("Configuring Jersey2 application");
    
    packages("com.bina.apac").register(JacksonFeature.class);
  }
}
