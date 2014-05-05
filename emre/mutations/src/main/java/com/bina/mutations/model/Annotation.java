package com.bina.mutations.model;


public class Annotation {
  
  public static enum Impact {
    DAMAGING, TOLERATED
  }

  public static enum Gene {
    BRCA1, BRCA2, KRAS, TPL1, PTEN
  }
  
  private Integer mutationId;
  private String transcriptId;
  private Gene gene;
  private Impact impact;

  public Annotation() {}

  public Annotation(Integer mutationId, String transcriptId, Gene gene, Impact impact) {
    this.mutationId = mutationId;
    this.transcriptId = transcriptId;
    this.gene = gene;
    this.impact = impact;
  }
  
  public Integer getMutationId() {
    return mutationId;
  }

  public void setMutationId(Integer mutationId) {
    this.mutationId = mutationId;
  }

  public String getTranscriptId() {
    return transcriptId;
  }

  public void setTranscriptId(String transcriptId) {
    this.transcriptId = transcriptId;
  }

  public Gene getGene() {
    return gene;
  }

  public void setGene(Gene gene) {
    this.gene = gene;
  }

  public Impact getImpact() {
    return impact;
  }

  public void setImpact(Impact impact) {
    this.impact = impact;
  }

}
