package com.bina.apac.resource;

import java.util.Map;

import javax.servlet.http.HttpServletRequest;
import javax.ws.rs.WebApplicationException;
import javax.ws.rs.client.Client;
import javax.ws.rs.client.ClientBuilder;
import javax.ws.rs.client.Entity;
import javax.ws.rs.client.Invocation;
import javax.ws.rs.core.Context;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.NewCookie;
import javax.ws.rs.core.Response;
import javax.ws.rs.core.Response.Status.Family;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.bina.apac.model.Login;
import com.fasterxml.jackson.databind.ObjectMapper;

public class BaseResource {

  private static final Logger logger = LoggerFactory.getLogger(BaseResource.class);
  
  protected static final String AP_REST_API_BASE_URL = "https://a.binacloud.com/api/v1";
  protected static final Client jerseyClient = ClientBuilder.newClient();
  protected static final ObjectMapper om = new ObjectMapper();
  
  @Context
  protected HttpServletRequest request;
 
  protected Response login() {
    return jerseyClient.target(AP_REST_API_BASE_URL)
                       .path("session")
                       .request(MediaType.APPLICATION_JSON)
                       .post(Entity.json(new Login("admin@bina.com", "BinaAdmin1")));
  }
  
  @SuppressWarnings("unchecked")
  protected void loginIfNecessary() {
    final Map<String, NewCookie> cookies = (Map<String, NewCookie>) request.getSession().getAttribute("cookies");
    if (cookies == null) {
      logger.info("Need to log in first. Loggin in...");
      
      final Response loginResponse = login();
      if (Family.SUCCESSFUL != loginResponse.getStatusInfo().getFamily()) {
        throw new WebApplicationException(Response.Status.UNAUTHORIZED);
      }
      
      request.getSession().setAttribute("cookies", loginResponse.getCookies());
    }
  }
 
  @SuppressWarnings("unchecked")
  protected Map<String,NewCookie> getCookiesFromSession(HttpServletRequest request) {
    return (Map<String,NewCookie>) request.getSession().getAttribute("cookies");
  }
  
  protected Invocation.Builder addCookiesToRequest(Invocation.Builder builder, Map<String, NewCookie> cookies) {
    for (final NewCookie nc : cookies.values()) {
      builder.cookie(nc.toCookie());
    }
    return builder;
  }
}
