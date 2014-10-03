package com.bina.apac.resource;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.zip.GZIPInputStream;

import javax.ws.rs.DefaultValue;
import javax.ws.rs.GET;
import javax.ws.rs.Path;
import javax.ws.rs.Produces;
import javax.ws.rs.QueryParam;
import javax.ws.rs.WebApplicationException;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;
import javax.ws.rs.core.Response.Status;

import org.apache.commons.lang.StringUtils;

import com.bina.apac.model.Activity;
import com.bina.apac.model.Notification;
import com.bina.apac.model.PagedResult;
import com.fasterxml.jackson.core.type.TypeReference;

@Produces(MediaType.APPLICATION_JSON)
@Path("activity")
public class ActivityResource extends BaseResource {

  @GET
  public List<Activity> getNotifications(@DefaultValue("0") @QueryParam("page") int page,
                                         @DefaultValue("50") @QueryParam("pagesize") int pageSize,
                                         @DefaultValue("DEBUG") @QueryParam("level") String level,
                                         @QueryParam("targetEmails") String targetEmails) {
    if (StringUtils.isBlank(targetEmails)) {
      throw new WebApplicationException(Response.status(Status.BAD_REQUEST)
                                        .entity("Please provide at least one email address")
                                        .build()); 
    }
    
    loginIfNecessary();

    final List<Activity> result = new ArrayList<>();
    try {
      String[] emailList = targetEmails.split("\\s*,\\s*");
      for (String email : emailList) {
        final Response response = getNotificationsForUser(email, level, page, pageSize); 
        final GZIPInputStream is = new GZIPInputStream(response.readEntity(InputStream.class));
        final PagedResult<Notification> pagedResult = om.readValue(is, new TypeReference<PagedResult<Notification>>() { });
        List<Notification> notifications = pagedResult != null ? pagedResult.getObjects() : null;
        if (notifications != null) {
          for (Notification n : notifications) {
            result.add(new Activity(email, n));
          } 
        }
      }
      
      Collections.sort(result, new Comparator<Activity>() {
        @Override
        public int compare(Activity first, Activity second) {
          return second.getDateCreated().compareTo(first.getDateCreated());
        }
      });
    } catch (Exception e) {
      throw new WebApplicationException(Response.status(Status.INTERNAL_SERVER_ERROR)
                                        .entity("Ooops... Something wrong happened")
                                        .build());
    }
    return (result.size() > pageSize ? result.subList(0, pageSize) : result);
  }
  
  private Response getNotificationsForUser(String email, String level, int page, int pageSize) {
    return addCookiesToRequest(jerseyClient.target(AP_REST_API_BASE_URL)
                                           .path("notification")
                                           .queryParam("page", page)
                                           .queryParam("pagesize", pageSize)
                                           .queryParam("targetEmail", email)
                                           .queryParam("level", level)
                                           .request(MediaType.APPLICATION_JSON),
                               getCookiesFromSession(request)).get();
  }
  
}
