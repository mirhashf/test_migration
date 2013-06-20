import java.io.Serializable;


public class Region implements Serializable {

  /**
   * 
   */
  private static final long serialVersionUID = -336449320692471875L;
  @Override
  public String toString() {
    return "Region [begins=" + begins + ", beginc=" + beginc + ", lasts=" + lasts + ", lastc="
        + lastc + ", numNormalReads=" + numNormalReads + "]";
  }
  public Region(String begins, Integer beginc, String lasts, Integer lastc, Integer numNormalReads) {
    super();
    this.begins = begins;
    this.beginc = beginc;
    this.lasts = lasts;
    this.lastc = lastc;
    this.numNormalReads = numNormalReads;
  }
  
  public String begins = "";
  public Integer beginc = -1;
  public String lasts = "";
  public Integer lastc = -1;
  public Integer numNormalReads=0;
  
  public boolean equals(Region r){
    return r.begins.equals(begins) && r.beginc==(beginc) && r.lasts.equals(lasts) && r.lastc==(lastc);
  }
  
  public int hashCode(){
    return begins.hashCode()+lasts.hashCode()+new Integer(beginc).hashCode()+new Integer(lastc).hashCode();
  }
  
}
