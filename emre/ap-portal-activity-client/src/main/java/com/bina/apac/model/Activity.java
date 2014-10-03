package com.bina.apac.model;

public class Activity {

  private String emailAddress;
  private Long id;
  private Long dateCreated;
  private String message;
  private String type;
  private String level;
  
  public Activity() {}
  
  public Activity(String emailAddress, Notification notification) {
    this.emailAddress = emailAddress;
    this.id = notification.getId();
    this.dateCreated = notification.getDateCreated();
    this.message = notification.getMessage();
    this.type = notification.getType();
    this.level = notification.getLevel();
  }
  
  public String getEmailAddress() {
    return emailAddress;
  }
  public void setEmailAddress(String emailAddress) {
    this.emailAddress = emailAddress;
  }
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
