import java.io.Serializable;
import java.util.HashMap;


public class CandidatorParameters implements Serializable {

  public String inputStats;
  public String inputBAM;
  public String scratchDir;
  public String chr;
  public Integer start;
  public Integer end;
  public int cut_sd;
  public int MINMAPQ;
  
  
  protected HashMap<String, Integer> chrLength;
  protected int numReads = 0;
  protected double chrDistanceToSeeAnomaly = Math.pow(10, 8);
  protected HashMap<String, ReadGroupInfo> readgroup;
  protected HashMap<String,Integer> numAllReadsPerLibrary=new HashMap<String, Integer>();
  protected HashMap<Integer,HashMap<String,Integer>> numAllBadReadsPerType=new HashMap<Integer, HashMap<String,Integer>>();
  protected HashMap<String,Double> readDensity=new HashMap<String,Double>();
  protected  double refLen;
  
  protected static final long serialVersionUID = 168412342202110909L;

}
