import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

import net.sf.samtools.SAMRecord;


public class RegionInfo implements Serializable {
  protected Integer reg_idx = 0;
  private static final long serialVersionUID = -9099909435879653037L;
  protected HashMap<Integer, ArrayList<SAMRecord>> regs = new HashMap<Integer, ArrayList<SAMRecord>>();
  protected HashMap<Integer, Region> regs_name = new HashMap<Integer, Region>();
  protected ArrayList<SAMRecord> reg_seq = new ArrayList<SAMRecord>();
  protected HashMap<String, ArrayList<Integer>> reads = new HashMap<String, ArrayList<Integer>>();
  
}
