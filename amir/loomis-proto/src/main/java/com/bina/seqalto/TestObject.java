// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import javax.persistence.Entity;

import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.Table;

/**
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 *
 */
@Entity
@Table(name = "TestObject")
public class TestObject {
  private int id;
  private String name;
  private int age;

  public void setId(int id) {
    this.id = id;
  }

  @Id
  @GeneratedValue
  public int getId() {
    return id;
  }

  public void setName(String name) {
    this.name = name;
  }

  public String getName() {
    return name;
  }

  public void setAge(int age) {
    this.age = age;
  }

  public int getAge() {
    return age;
  }
}
