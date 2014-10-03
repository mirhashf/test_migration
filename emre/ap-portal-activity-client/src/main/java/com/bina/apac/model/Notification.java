package com.bina.apac.model;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;

@JsonIgnoreProperties(ignoreUnknown = true)
public class Notification {

  private Long id;
  private Long dateCreated;
  private String message;
  private String type;
  private String level;
  
  public Long getId() {
    return id;
  }
  public void setId(Long id) {
    this.id = id;
  }
  public Long getDateCreated() {
    return dateCreated;
  }
  public void setDateCreated(Long dateCreated) {
    this.dateCreated = dateCreated;
  }
  public String getMessage() {
    return message;
  }
  public void setMessage(String message) {
    this.message = message;
  }
  public String getType() {
    return type;
  }
  public void setType(String type) {
    this.type = type;
  }
  public String getLevel() {
    return level;
  }
  public void setLevel(String level) {
    this.level = level;
  }
  
}
