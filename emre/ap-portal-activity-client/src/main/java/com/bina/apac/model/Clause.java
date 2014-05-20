package com.bina.apac.model;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;

@JsonIgnoreProperties(ignoreUnknown = true)
public class Clause {

  private String dataSourceName;
  private String columnName;
  private String operator;
  private String operand;
  
  public String getDataSourceName() {
    return dataSourceName;
  }
  public void setDataSourceName(String dataSourceName) {
    this.dataSourceName = dataSourceName;
  }
  public String getColumnName() {
    return columnName;
  }
  public void setColumnName(String columnName) {
    this.columnName = columnName;
  }
  public String getOperator() {
    return operator;
  }
  public void setOperator(String operator) {
    this.operator = operator;
  }
  public String getOperand() {
    return operand;
  }
  public void setOperand(String operand) {
    this.operand = operand;
  }
  
}
