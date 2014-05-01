package com.bina.apac.model;

import java.util.ArrayList;
import java.util.List;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;

@JsonIgnoreProperties(ignoreUnknown = true)
public class Query {

  private String dataSetName;
  private String sampleName;
  private List<Clause> andClauses = new ArrayList<>();
  
  public String getDataSetName() {
    return dataSetName;
  }
  public void setDataSetName(String dataSetName) {
    this.dataSetName = dataSetName;
  }
  public String getSampleName() {
    return sampleName;
  }
  public void setSampleName(String sampleName) {
    this.sampleName = sampleName;
  }
  public List<Clause> getAndClauses() {
    return andClauses;
  }
  public void setAndClauses(List<Clause> andClauses) {
    this.andClauses = andClauses;
  }
  
}
