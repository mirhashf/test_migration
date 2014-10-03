package com.bina.apac.resource;

import javax.ws.rs.GET;
import javax.ws.rs.Path;
import javax.ws.rs.PathParam;
import javax.ws.rs.Produces;
import javax.ws.rs.core.MediaType;

import com.bina.apac.model.Query;

@Produces(MediaType.APPLICATION_JSON)
@Path("queries")
public class QueryResource extends BaseResource {

  @GET
  @Path("{id}")
  public Query getQueryById(@PathParam("id") long id) {
    loginIfNecessary();
    
    return addCookiesToRequest(jerseyClient.target(AP_REST_API_BASE_URL)
                               .path(String.format("query/%d", id))
                               .request(MediaType.APPLICATION_JSON),
                               getCookiesFromSession(request)).get(Query.class);
  }
}
