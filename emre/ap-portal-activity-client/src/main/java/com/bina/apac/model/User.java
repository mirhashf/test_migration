package com.bina.apac.model;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;

@JsonIgnoreProperties(ignoreUnknown = true)
public class User {

  private Long id;
  private String emailAddress;
  private String firstName;
  private String lastName;
  private String organization;
  private boolean emailVerified;
  private Long dateCreated;
  
  public Long getId() {
    return id;
  }
  public void setId(Long id) {
    this.id = id;
  }
  public String getEmailAddress() {
    return emailAddress;
  }
  public void setEmailAddress(String emailAddress) {
    this.emailAddress = emailAddress;
  }
  public String getFirstName() {
    return firstName;
  }
  public void setFirstName(String firstName) {
    this.firstName = firstName;
  }
  public String getLastName() {
    return lastName;
  }
  public void setLastName(String lastName) {
    this.lastName = lastName;
  }
  public String getOrganization() {
    return organization;
  }
  public void setOrganization(String organization) {
    this.organization = organization;
  }
  public boolean isEmailVerified() {
    return emailVerified;
  }
  public void setEmailVerified(boolean emailVerified) {
    this.emailVerified = emailVerified;
  }
  public Long getDateCreated() {
    return dateCreated;
  }
  public void setDateCreated(Long dateCreated) {
    this.dateCreated = dateCreated;
  }
 
}
