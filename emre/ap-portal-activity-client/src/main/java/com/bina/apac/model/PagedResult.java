package com.bina.apac.model;

import java.util.List;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;

@JsonIgnoreProperties(ignoreUnknown = true)
public class PagedResult<T> {

  private List<T> objects;
  private int page;
  private int pageSize;
  private int totalPages;
  private long total;
  
  public List<T> getObjects() {
    return objects;
  }
  public void setObjects(List<T> objects) {
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
  public int getTotalPages() {
    return totalPages;
  }
  public void setTotalPages(int totalPages) {
    this.totalPages = totalPages;
  }
  public long getTotal() {
    return total;
  }
  public void setTotal(long total) {
    this.total = total;
  }
  
}
