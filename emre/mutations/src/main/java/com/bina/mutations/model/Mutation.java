package com.bina.mutations.model;


public class Mutation {

  public static enum Base {
    A, C, T, G
  }

  public static enum Chromosome {
    CHR1, CHR2, CHR3, CHR4, CHR5, CHR6, CHR7, CHR8, CHR9, CHR10, CHR11, CHR12, CHR13, 
    CHR14, CHR15, CHR16, CHR17, CHR18, CHR19, CHR20, CHR21, CHRX, CHRY, CHRM
  }

  private Integer id;
  private Chromosome chromosome;
  private Integer position;
  private String reference;
  private String alternate;
  private Integer readDepth;

  public Mutation() {}

  public Mutation(Integer id,
                  Chromosome chromosome,
                  Integer position,
                  String reference,
                  String alternate,
                  Integer readDepth) {
    this.id = id;
    this.chromosome = chromosome;
    this.position = position;
    this.reference = reference;
    this.alternate = alternate;
    this.readDepth = readDepth;
  }

  public Integer getId() {
    return id;
  }

  public void setId(Integer id) {
    this.id = id;
  }

  public Chromosome getChromosome() {
    return chromosome;
  }

  public void setChromosome(Chromosome chromosome) {
    this.chromosome = chromosome;
  }

  public Integer getPosition() {
    return position;
  }

  public void setPosition(Integer position) {
    this.position = position;
  }

  public String getReference() {
    return reference;
  }

  public void setReference(String reference) {
    this.reference = reference;
  }

  public String getAlternate() {
    return alternate;
  }

  public void setAlternate(String alternate) {
    this.alternate = alternate;
  }

  public Integer getReadDepth() {
    return readDepth;
  }

  public void setReadDepth(Integer readDepth) {
    this.readDepth = readDepth;
  }

}
