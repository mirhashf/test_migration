package com.bina.mutations.model;

import java.util.ArrayList;
import java.util.List;

public class QueryResult {

  private List<?> objects = new ArrayList<>();
  private int page;
  private int pageSize;
  private int total;
 
  public QueryResult() {}
  
  public QueryResult(List<?> objects, int page, int pageSize, int total) {
    this.objects = objects;
    this.page = page;
    this.pageSize = pageSize;
    this.total = total;
  }

  public List<?> getObjects() {
    return objects;
  }
  public void setObjects(List<?> objects) {
    this.objects = objects;
  }
  public int getPage() {
    return page;
  }
  public void setPage(int page) {
    this.page = page;
  }
  public int getPageSize() {
    return pageSize;
  }
  public void setPageSize(int pageSize) {
    this.pageSize = pageSize;
  }
  public int getTotal() {
    return total;
  }
  public void setTotal(int total) {
    this.total = total;
  }
   
}