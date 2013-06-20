import java.io.Serializable;

import net.sf.samtools.SAMRecord;

public class ReadGroupInfo implements Serializable {
  /**
   * 
   */
  private static final long serialVersionUID = 7692367297766616475L;
  public double mean;

  @Override
  public String toString() {
    return "readGroupInfo [mean=" + mean + ", std=" + std + ", readlen=" + readlen + ", count="
        + count + "]";
  }

  public double std;
  public double readlen;
  public double count;
  private int MINMAPQ;

  public double upper;
  public double lower;
  
  public ReadGroupInfo(int MINMAPQ) {
    this.MINMAPQ = MINMAPQ;
  }

  public void addInfo(SAMRecord record) {

    int insertSize = Math.abs(record.getInferredInsertSize());
    if (record.getMappingQuality() < MINMAPQ || record.getDuplicateReadFlag()
        || record.getMateUnmappedFlag()) return;
    if (insertSize > 2000) System.out.println(insertSize);
    mean += insertSize;
    std += insertSize * insertSize;
    readlen += record.getReadLength();
    count++;
  }

}
