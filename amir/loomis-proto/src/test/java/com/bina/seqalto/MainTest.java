package com.bina.seqalto;


import org.codehaus.jackson.jaxrs.JacksonJsonProvider;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;
import org.springframework.web.context.ContextLoaderListener;
import org.springframework.web.context.request.RequestContextListener;

import com.sun.jersey.api.client.WebResource;
import com.sun.jersey.api.client.config.ClientConfig;
import com.sun.jersey.api.client.config.DefaultClientConfig;
import com.sun.jersey.api.json.JSONConfiguration;
import com.sun.jersey.spi.spring.container.servlet.SpringServlet;
import com.sun.jersey.test.framework.JerseyTest;
import com.sun.jersey.test.framework.WebAppDescriptor;


@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration(locations = {"classpath:applicationContext.xml"})
public class MainTest extends JerseyTest {

  private static final Logger logger = LoggerFactory.getLogger(MainTest.class);

  WebResource resource;


  public MainTest() {
    // took me a while to figure these out...
    // i wanted to make it work with in-memory container but seems like this is
    // the only option
    super(new WebAppDescriptor.Builder("com.bina.seqalto")
        .contextParam("contextConfigLocation", "classpath:applicationContext.xml")
        .contextPath("/root").servletClass(SpringServlet.class)
        .contextListenerClass(ContextLoaderListener.class)
        .requestListenerClass(RequestContextListener.class).clientConfig(createClientConfig())
        .initParam(JSONConfiguration.FEATURE_POJO_MAPPING, "true").build());
  }

  private static ClientConfig createClientConfig() {
    ClientConfig config = new DefaultClientConfig();
    config.getClasses().add(JacksonJsonProvider.class);
    config.getFeatures().put(JSONConfiguration.FEATURE_POJO_MAPPING, true);
    return config;
  }
  
 
  @Test
  public void testZooKeeper(){
    logger.info("Testing ZooKeeper");
  }
}
