package tmp;
import java.io.Serializable;


class Coveragem implements Serializable {


  private static final long serialVersionUID = 1661574348164780923L;

  public Coveragem() {

  }

  public Coveragem(long numReads, double readDensity, double total_doc, double insertsize_doc,
      double readLen) {
    super();
    this.numReads = numReads;
    this.readDensity = readDensity;
    this.total_doc = total_doc;
    this.insertsize_doc = insertsize_doc;
    this.readLen = readLen;
  }

  @Override
  public String toString() {
    return "Coverage [numReads=" + numReads + ", readDensity=" + readDensity + ", readLen="
        + readLen + ", eventSize=" + eventSize + ", total_doc=" + total_doc + ", insertsize_doc="
        + insertsize_doc + "]";
  }

  public double readLen;
  public long numReads;
  public double insertsize;

  public double readDensity;
  public double total_doc;
  public double insertsize_doc;
  public double eventSize;
}
