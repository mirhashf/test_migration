// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import javax.ws.rs.WebApplicationException;
import javax.ws.rs.core.Response;

import com.sun.jersey.spi.container.ContainerRequest;
import com.sun.jersey.spi.container.ContainerRequestFilter;

/**
 * Intercepter for all REST calls to make sure it is coming to a Leader. Otherwise it returns a
 * Forbidden status.
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com) Created on Feb 9, 2012
 */
public class JerseyInterceptor implements ContainerRequestFilter {
  /*
   * (non-Javadoc)
   * 
   * @see com.sun.jersey.spi.container.ContainerRequestFilter#filter(com.sun.jersey.spi.container.
   * ContainerRequest)
   */
  public ContainerRequest filter(ContainerRequest request) {
    if (ZooKeeperService.getIsLeader()) {
      return request;
    } else {
      throw new WebApplicationException(Response.Status.FORBIDDEN);
    }
  }
}
