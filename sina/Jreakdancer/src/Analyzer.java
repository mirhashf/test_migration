import java.io.*;
import java.util.*;

import net.sf.samtools.BAMIndexer;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;


public class Analyzer {

  private String scratchDir;

  protected AnalyzerParameters ap;
  protected CandidatorParameters cp;

  public Analyzer(String scratchDir, int minReadPair, int scoreThreshold) {
    super();
    this.scratchDir = scratchDir;
    ap = new AnalyzerParameters(minReadPair, scoreThreshold);
    try {
      getCandidatorParameter();
    } catch (Exception e) {
      e.printStackTrace();
    }
    File scratchDirFile = new File(cp.scratchDir);
    if (!scratchDirFile.exists()) scratchDirFile.mkdirs();
    analysis();
    System.out.println(ap.numEvents+" events detected");
  }

  private void analysis() {

    SAMFileReader reader = new SAMFileReader(new File(cp.inputBAM));
    SAMRecordIterator iterator = reader.iterator();
    if (cp.chr.trim().length() > 0) {
      if (!reader.hasIndex())
        BAMIndexer.createAndWriteIndex(new File(cp.inputBAM), new File(cp.inputBAM + "bai"), false);
      iterator.close();
      iterator = reader.query(cp.chr, (int) cp.start, (int) cp.end, true);
    }

    while (iterator.hasNext()) {
      if (ap.numReads > Candidator.numReadsToProcess) break;
      SAMRecord record = iterator.next();
      ap.numReads++;
      analyze(record);
    }
    if (ap.regionInfo.reg_seq.size() != 0) badRegionProcess(false);
    buildConnection();
    iterator.close();
    reader.close();
  }


  private void buildConnection() {
    RegionGraph graph = new RegionGraph();
    for (String read_name : ap.regionInfo.reads.keySet()) {
      ArrayList<Integer> regionPairList = ap.regionInfo.reads.get(read_name);
      if (regionPairList.size() != 2) {
        /*System.err
            .println("the readname does not correspond to a pair of regions (maybe it is not paired)"+regionPairList.size());
       */
       // System.out.println(read_name+"\t"+ap.regionInfo.regs.get(0)+"\t");
        //System.exit(1);
        continue;
      }
      graph.add(regionPairList.get(0), regionPairList.get(1));
    }

    HashSet<Integer> freeRegions = new HashSet();

    while(graph.keySet().size()>0) {
      Integer node=graph.keySet().iterator().next();
      LinkedList<Integer> bfsNodes = new LinkedList<Integer>();
      bfsNodes.add(node);
      while (bfsNodes.size() > 0) {
        LinkedList<Integer> newBfsNodes = new LinkedList<Integer>();
        for (Integer node0 : bfsNodes) {
          if(!graph.containsKey(node0))
            continue;
          Iterator<Integer> iter=graph.get(node0).keySet().iterator();
          while(iter.hasNext()){
            Integer node1=iter.next();
            int numLinks = graph.get(node0).get(node1);
            if (numLinks < ap.minReadPair) continue;
            SVDetector detector = new SVDetector(this, freeRegions);
            detector.examinePairOfRegions(node0, node1);
            graph.remove(node0, node1);
            newBfsNodes.add(node1);
          }
          
          graph.remove(node0);
        }
        bfsNodes = newBfsNodes;
      }
    }

    ap.freeRegions(freeRegions);
  }

 
  private void badRegionProcess(boolean inMiddle) {
    double seq_coverage =
        (double) ap.ntotal_nucleotids / ap.maxReadLength;
    if (ap.lastc - ap.beginc > ap.minEventSize && seq_coverage < ap.SEQ_COVERAGE_LIM)
      startLargeBadRegionDetection(inMiddle);
    else
      startSmallBadRegionDetection(inMiddle);

  }

  private void analyze(SAMRecord record) {
    double flag = AlnParser.alnParse(record);
    String rg = record.getReadGroup().getId();
    if (record.getReadLength() > ap.maxReadLength) ap.maxReadLength = record.getReadLength();
    if (record.getMappingQuality() > cp.MINMAPQ && flag < 32 && flag >= 18) update_readcounts(rg);
    if (record.getMappingQuality() < cp.MINMAPQ || flag == 0 || record.getDuplicateReadFlag())
      return;

    int insertSize = record.getInferredInsertSize();
    flag = updateFlag(flag, rg, insertSize);

    if (flag < 32 && Math.abs(insertSize) > ap.maxEventSize) return;
    if (flag == 18) {
      if (ap.normal_switch == true && insertSize > 0) ap.nnormalReads++;
      return;
    }
    record.setAttribute("fz", (int)flag);
    
    if (ap.normal_switch == true) ap.ntotal_nucleotids += record.getReadLength();

    boolean isNewBadRegion = false;
    if (!record.getReferenceName().equals(ap.lasts)
        || record.getAlignmentStart() - ap.lastc > cp.chrDistanceToSeeAnomaly)
      isNewBadRegion = true;

    if (isNewBadRegion == true) {
      badRegionProcess(true);

      ap.endOfBadRegionResetParamters();
      ap.beginc = record.getAlignmentStart();
      ap.begins = record.getReferenceName();

    }
    ap.endOfBadReadResetParameters(record);

  }

  private void startLargeBadRegionDetection(boolean inMiddle) {
    Region newRegion = new Region(ap.begins, ap.beginc, ap.lasts, ap.lastc, ap.nnormalReads);
    Integer newRegionIndex = ap.regionInfo.reg_idx;
    ap.regionInfo.reg_idx++;
    ap.regionInfo.regs_name.put(newRegionIndex, newRegion);

    if (inMiddle == true) updateRegionReadCount(newRegionIndex);

    ArrayList<SAMRecord> regionRecords = new ArrayList<SAMRecord>();
    for (SAMRecord record : ap.regionInfo.reg_seq) {
      regionRecords.add(record);
      if (!ap.regionInfo.reads.containsKey(record.getReadName()))
        ap.regionInfo.reads.put(record.getReadName(), new ArrayList<Integer>());
      ap.regionInfo.reads.get(record.getReadName()).add(newRegionIndex);
    }
    ap.regionInfo.regs.put(newRegionIndex, regionRecords);

    ap.idx_buff++;
    if (ap.idx_buff > ap.BUFFER_SIZE) {
      buildConnection();
      ap.idx_buff = 0;
    }

  }

  private void updateRegionReadCount(Integer newRegionIndex) {
    if (!ap.rcInfo.read_count_roi_map.containsKey(newRegionIndex))
      ap.rcInfo.read_count_roi_map.put(newRegionIndex, new HashMap<String, Integer>());
    for (String rg : ap.rcInfo.nread_roi.keySet()) {
      if (!ap.rcInfo.read_count_roi_map.get(newRegionIndex).containsKey(rg))
        ap.rcInfo.read_count_roi_map.get(newRegionIndex).put(rg, (Integer) 0);
      ap.rcInfo.read_count_roi_map.get(newRegionIndex).put(rg,
          ap.rcInfo.read_count_roi_map.get(newRegionIndex).get(rg) + ap.rcInfo.nread_roi.get(rg));
    }
    if (!ap.rcInfo.read_cound_fr_map.containsKey(newRegionIndex))
      ap.rcInfo.read_cound_fr_map.put(newRegionIndex, new HashMap<String, Integer>());
    for (String rg : ap.rcInfo.nread_fr.keySet()) {
      if (!ap.rcInfo.read_cound_fr_map.get(newRegionIndex).containsKey(rg))
        ap.rcInfo.read_cound_fr_map.get(newRegionIndex).put(rg, (Integer) 0);
      Integer diff = ap.rcInfo.nread_fr.get(rg);
      if (ap.rcInfo.read_count_roi_map.containsKey(newRegionIndex)
          && ap.rcInfo.read_count_roi_map.get(newRegionIndex).containsKey(rg))
        diff -= ap.rcInfo.read_count_roi_map.get(newRegionIndex).get(rg);
      ap.rcInfo.read_cound_fr_map.get(newRegionIndex).put(rg, diff);
    }

  }

  private void startSmallBadRegionDetection(boolean inMiddle) {
    if (inMiddle == true) {
      if (ap.regionInfo.reg_idx >= 1) {
        if (!ap.rcInfo.read_count_roi_map.containsKey(ap.regionInfo.reg_idx - 1))
          ap.rcInfo.read_count_roi_map.put(ap.regionInfo.reg_idx - 1,
              new HashMap<String, Integer>());
        for (String rg : ap.rcInfo.possible_fake_data.keySet()) {
          if (!ap.rcInfo.read_count_roi_map.get(ap.regionInfo.reg_idx - 1).containsKey(rg))
            ap.rcInfo.read_count_roi_map.get(ap.regionInfo.reg_idx - 1).put(rg, (Integer) 0);
          ap.rcInfo.read_count_roi_map.get(ap.regionInfo.reg_idx - 1).put(
              rg,
              ap.rcInfo.read_count_roi_map.get(ap.regionInfo.reg_idx - 1).get(rg)
                  + ap.rcInfo.possible_fake_data.get(rg));
        }
      }
    }

    if (ap.regionInfo.reg_seq.size() > 0) for (SAMRecord record : ap.regionInfo.reg_seq)
      ap.regionInfo.reads.remove(record.getReadName());


  }

  private double updateFlag(double flag, String rg, int insertSize) {
    insertSize = Math.abs(insertSize);
    ReadGroupInfo info = cp.readgroup.get(rg);
    if (insertSize > info.upper && flag == 18) flag = 2;
    if (insertSize <= info.upper && flag == 2) flag = 18;

    if (insertSize < info.lower && flag == 18) flag = 3;
    if (insertSize >= info.lower && flag == 3) flag = 18;
    if (flag == 20) flag = 4;
    if (flag == 8) flag = 1;
    return flag;
  }

  private void update_readcounts(String rg) {
    if (!ap.rcInfo.nread_roi.containsKey(rg)) ap.rcInfo.nread_roi.put(rg, (Integer) 0);
    ap.rcInfo.nread_roi.put(rg, ap.rcInfo.nread_roi.get(rg) + 1);

    if (!ap.rcInfo.nread_fr.containsKey(rg)) ap.rcInfo.nread_fr.put(rg, (Integer) 0);
    ap.rcInfo.nread_fr.put(rg, ap.rcInfo.nread_fr.get(rg) + 1);

    if (!ap.rcInfo.possible_fake_data.containsKey(rg))
      ap.rcInfo.possible_fake_data.put(rg, (Integer) 0);
    ap.rcInfo.possible_fake_data.put(rg, ap.rcInfo.possible_fake_data.get(rg) + 1);

  }



  private void getCandidatorParameter() throws Exception {
    ObjectInputStream r =
        new ObjectInputStream(new FileInputStream(new File(scratchDir,
            Candidator.CANDIDATE_PARAMETER_FILE).getAbsolutePath()));
    cp = (CandidatorParameters) r.readObject();
    r.close();

  }

  public static void getAnalyzer(String[] args) throws IOException {
    int count = 0;
    String inputBam = args[count++];
    String inputStats = args[count++];
    String scratchDir = args[count++];
    String chr = "";
    if (args.length > count) chr = args[count++];
    Integer start = 0;
    if (args.length > count) start = Integer.valueOf(args[count++]);
    Integer end = 0;
    if (args.length > count) start = Integer.valueOf(args[count++]);
    int cut_sd = 3;
    if (args.length > count) cut_sd = Integer.valueOf(args[count++]);
    int minMapQ = 35;
    if (args.length > count) minMapQ = Integer.valueOf(args[count++]);
    int minReadPair = 2;
    if (args.length > count) minReadPair = Integer.valueOf(args[count++]);
    int scoreThreshold = 29;
    if (args.length > count) scoreThreshold = Integer.valueOf(args[count++]);
    Analyzer candidator = new Analyzer(scratchDir, minReadPair, scoreThreshold);

  }
}
