package com.bina.mutations.filter;

import java.io.IOException;

import javax.servlet.Filter;
import javax.servlet.FilterChain;
import javax.servlet.FilterConfig;
import javax.servlet.ServletException;
import javax.servlet.ServletRequest;
import javax.servlet.ServletResponse;
import javax.servlet.http.HttpServletResponse;

public class CORSFilter implements Filter {

  private static final String ACCESS_CONTROL_ALLOW_ORIGIN = "Access-Control-Allow-Origin";
  private static final String ACCESS_CONTROL_ALLOW_METHODS = "Access-Control-Allow-Methods";
  private static final String ALLOWED_METHODS = "POST, GET, OPTIONS, DELETE, PUT";
  private static final String ACCESS_CONTROL_ALLOW_HEADERS = "Access-Control-Allow-Headers";
  private static final String ALLOWED_HEADERS = "X-Requested-With, Content-Type";

  @Override
  public void init(FilterConfig config) throws ServletException {}
  
  @Override
  public void destroy() {}

  @Override
  public void doFilter(ServletRequest request, ServletResponse response,
                       FilterChain filterChain) throws IOException, ServletException {
      HttpServletResponse resp = (HttpServletResponse) response;

      resp.addHeader(ACCESS_CONTROL_ALLOW_ORIGIN, "*");
      resp.addHeader(ACCESS_CONTROL_ALLOW_METHODS, ALLOWED_METHODS);
      resp.addHeader(ACCESS_CONTROL_ALLOW_HEADERS, ALLOWED_HEADERS);

      filterChain.doFilter(request, response);
  }

}