import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import net.sf.samtools.SAMRecord;

public class AnalyzerParameters implements Serializable {

  public static int numEvents; 
  private static final long serialVersionUID = -1912921259452264411L;

  protected ReadCountInfo rcInfo = new ReadCountInfo();
  protected RegionInfo regionInfo = new RegionInfo();

  protected String begins = "";
  protected Integer beginc = -1;
  protected String lasts = "";
  protected Integer lastc = -1;

  protected int numReads = 0;
  protected int nnormalReads = 0;
  protected Integer ntotal_nucleotids = 0; // global

  protected int minEventSize = 100; // minimum region length
  protected double maxEventSize = Math.pow(10, 9);
  public int maxReadLength = 100;

  protected int scoreThreshold;
  protected int minReadPair;

  protected boolean normal_switch = false;
  protected final int SEQ_COVERAGE_LIM = 1000;
  protected final int BUFFER_SIZE = 100;
  protected int idx_buff = 0;


  public AnalyzerParameters(int minReadPair, int scoreThreshold) {
    super();
    this.minReadPair = minReadPair;
    this.scoreThreshold = scoreThreshold;
  }

  protected static HashMap<Integer, String> SVtype;
  static {
    SVtype = new HashMap<Integer, String>();
    SVtype.put(1, "INV");
    SVtype.put(2, "DEL");
    SVtype.put(3, "INS");
    SVtype.put(4, "ITX");
    SVtype.put(8, "INV");
    SVtype.put(32, "CTX");
  }

  public void endOfBadRegionResetParamters() {
    regionInfo.reg_seq.clear();
    normal_switch = false;
    nnormalReads = 0;
    ntotal_nucleotids = 0;

    rcInfo.possible_fake_data.clear();
    rcInfo.nread_roi.clear();
    rcInfo.nread_fr.clear();
  }

  public void endOfBadReadResetParameters(SAMRecord record) {

    regionInfo.reg_seq.add(record);
    if (regionInfo.reg_seq.size() == 1) normal_switch = true;
    lasts = record.getReferenceName();
    lastc = record.getAlignmentStart();
    rcInfo.nread_roi.clear();

  }

  public void freeRegions(HashSet<Integer> freeRegions) {
    for (Integer freeRegion : freeRegions) {
      ArrayList<SAMRecord> freeReads = regionInfo.regs.get(freeRegion);
      if (freeReads.size() < minReadPair) {
        for (SAMRecord record : freeReads)
          regionInfo.reads.remove(record.getReadName());
        regionInfo.regs.remove(freeRegion);
        regionInfo.regs_name.remove(freeRegion);
      }
    }
  }
}


class ReadCountInfo implements Serializable {

  private static final long serialVersionUID = -8950452201267234089L;
  protected HashMap<String, Integer> nread_roi = new HashMap<String, Integer>();
  protected HashMap<String, Integer> nread_fr = new HashMap<String, Integer>();
  protected HashMap<String, Integer> possible_fake_data = new HashMap<String, Integer>();

  protected HashMap<Integer, HashMap<String, Integer>> read_count_roi_map =
      new HashMap<Integer, HashMap<String, Integer>>();
  protected HashMap<Integer, HashMap<String, Integer>> read_cound_fr_map =
      new HashMap<Integer, HashMap<String, Integer>>();
}
