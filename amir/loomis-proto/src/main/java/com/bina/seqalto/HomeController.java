// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import javax.annotation.Resource;
import javax.inject.Inject;
import javax.ws.rs.GET;
import javax.ws.rs.Path;
import javax.ws.rs.PathParam;
import javax.ws.rs.Produces;
import javax.ws.rs.core.MediaType;

import org.hibernate.Session;
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


  @GET
  @Path("new_job/{name}")
  @Produces(MediaType.APPLICATION_JSON)
  public TestObject setIt(@PathParam("name") String name) {
    Session session = getSessionFactory().openSession();
    session.beginTransaction();
    TestObject object = new TestObject();
    object.setName(name);
    session.save(object);
    session.getTransaction().commit();
    return object;
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
