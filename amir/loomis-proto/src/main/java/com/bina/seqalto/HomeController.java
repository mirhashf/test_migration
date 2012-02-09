// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import javax.annotation.Resource;
import javax.inject.Inject;
import javax.ws.rs.GET;
import javax.ws.rs.Path;
import javax.ws.rs.PathParam;
import javax.ws.rs.Produces;
import javax.ws.rs.core.MediaType;

import org.hibernate.SessionFactory;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

/**
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 * 
 * Stateless bean to respond to JSON API Calls.
 * Built on top of Jersey
 */
@Path("/")
@Component
@Scope("request")
public class HomeController {
  @Resource(name = "sessionFactory")
  @Inject
  public SessionFactory sessionFactory;
  @Inject
  public Leader leader;

  @GET
  @Path("new_job/{binary}")
  @Produces(MediaType.APPLICATION_JSON)
  public Object newJob(@PathParam("binary") String binary) throws Exception {
    if(ZooKeeperService.getIsLeader()){
      JobRequest job = new JobRequest();
      job.setBinary(binary);
      leader.process(job);
      return job;
    }else{
      return "I'm not a leader! Please contact leader...";
    }
  }

  @GET
  @Path("goodbye")
  @Produces(MediaType.APPLICATION_JSON)
  public TestObject fuckIt() {
    TestObject test = new TestObject();
    return test;
  }

  @GET
  @Path("areyouleader")
  @Produces(MediaType.APPLICATION_JSON)
  public boolean areYouLeader() {
    return ZooKeeperService.getIsLeader();
  }


  public void setSessionFactory(SessionFactory sessionFactory) {
    this.sessionFactory = sessionFactory;
  }

  public SessionFactory getSessionFactory() {
    return sessionFactory;
  }
}
