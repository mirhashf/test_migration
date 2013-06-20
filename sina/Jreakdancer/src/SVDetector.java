import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import net.sf.samtools.SAMRecord;


public class SVDetector {

  private int flag;
  private SVBreakPointParameters svbp;

  private int nread_pairs = 0;
  private HashMap<String, SAMRecord> readPair;
  private HashMap<Integer, Integer> SVTypeReadcount; // first one is flags
  protected HashMap<Integer, HashMap<String, Integer>> SVTypeLibraryReadcount;
  private ArrayList<HashMap<String, Integer>> SVTypeOrientCounts;
  private HashMap<Integer, HashMap<String, Integer>> SVTypeLibraryMeanspan;// diff span distance;
  private ArrayList<SAMRecord> supportReads;
  private ArrayList<String> freeReads = new ArrayList<String>();

  public Analyzer an;
  private HashSet<Integer> freeRegions;

  public SVDetector(Analyzer an, HashSet<Integer> freeRegions) {
    this.an = an;
    this.freeRegions = freeRegions;
    readPair = new HashMap<String, SAMRecord>();
    SVTypeReadcount = new HashMap<Integer, Integer>();
    SVTypeLibraryReadcount = new HashMap<Integer, HashMap<String, Integer>>();
    SVTypeOrientCounts = new ArrayList<HashMap<String, Integer>>();
    SVTypeLibraryMeanspan = new HashMap<Integer, HashMap<String, Integer>>();
    supportReads = new ArrayList<SAMRecord>();
  }

  public void examinePairOfRegions(Integer node0, Integer node1) {
    examinRegion(node0);
    examinRegion(node1);

    cleanNonSupportives(node0);
    cleanNonSupportives(node1);

    if (nread_pairs >= an.ap.minReadPair) {
      flag = calculateMaxScore();
      if (SVTypeReadcount.get(flag) >= an.ap.minReadPair) {
        getEventBreakPoints(node0, node1);
        removeFreeReads(node0, node1);
      }
    }

  }

  private void removeFreeReads(Integer node0, Integer node1) {
    for (String record : freeReads)
      an.ap.regionInfo.reads.remove(record);
    freeRegions.add(node0);
    freeRegions.add(node1);

  }

  private void getEventBreakPoints(Integer node0, Integer node1) {
    svbp = new SVBreakPointParameters();
    getEventBreakPoints(node0);
    getEventBreakPoints(node1);
    getCNV(node0, node1);
    getScore(node0,node1);
    if (svbp.PhredQ > an.ap.scoreThreshold) {
      an.ap.numEvents++;
      printSV();
    }
  }


  private void printSV() {
    System.out.print(svbp.sv_chr1 + "\t" + svbp.sv_pos1 + "\t" + svbp.sv_ori1 + "\t");
    System.out.print(svbp.sv_chr2 + "\t" + svbp.sv_pos2 + "\t" + svbp.sv_ori2 + "\t");
    System.out.print(svbp.SVType + "\t" + (svbp.sv_pos2-svbp.sv_pos1)+"\t"+ svbp.PhredQ + "\t" + SVTypeReadcount.get(flag));
    System.out.println();

  }

  private void getScore(Integer node0,Integer node1) {
    ArrayList<Integer> regions=new ArrayList<Integer>();
    regions.add(node0); regions.add(node1);
    ScoreEstimator estimator=new ScoreEstimator(this);
    double LogPvalue = estimator.computeProbScore(regions,flag);
    svbp.PhredQ = (int) Math.round(-10 * LogPvalue / Math.log(10));
    if (svbp.PhredQ > 99) svbp.PhredQ = 99;
    if (an.ap.SVtype.containsKey(flag)) svbp.SVType = an.ap.SVtype.get(flag);

    // make the coordinates with base 1
    svbp.sv_pos1 = svbp.sv_pos1 + 1;
    svbp.sv_pos2 = svbp.sv_pos2 + 1;

  }

  private void getCNV(Integer node0,Integer node1) {
    HashMap<String, Double> copy_number = new HashMap<String, Double>();
   
      for (String rg : svbp.read_count.keySet()) {
      copy_number.put(rg,
          2 * (double) svbp.read_count.get(rg)
              / ((double) an.cp.readDensity.get(rg) * (svbp.sv_pos2 - svbp.sv_pos1)));
      svbp.copy_number += copy_number.get(rg);
    }
    if(svbp.read_count.size()>0)
    svbp.copy_number /= (2.0 * (double) svbp.read_count.size());

    if (flag != 4 && flag != 8 && svbp.sv_pos1 + an.ap.maxReadLength - 5 < svbp.sv_pos2)
      svbp.sv_pos1 += an.ap.maxReadLength - 5;

  }

  private void getEventBreakPoints(Integer node) {

    String chr = an.ap.regionInfo.regs_name.get(node).begins;
    Integer start = an.ap.regionInfo.regs_name.get(node).beginc;
    Integer end = an.ap.regionInfo.regs_name.get(node).lastc;
    Integer nrp = an.ap.regionInfo.regs_name.get(node).numNormalReads;

    if (SVTypeOrientCounts.size() == 0) {
      System.out.println("SVTypeOrientCounts is empty! going to halt!");
      System.exit(1);
    }
    HashMap<String, Integer> ori_readcount = SVTypeOrientCounts.get(0);
    SVTypeOrientCounts.remove(0);


    if (!svbp.sv_chr1.equals("") && !svbp.sv_chr2.equals("")) {
      switch (flag) {
        case 4:
          svbp.sv_pos2 = end + an.ap.maxReadLength - 5;
          break;
        case 1:
          svbp.sv_pos1 = svbp.sv_pos2;
          svbp.sv_pos2 = end + an.ap.maxReadLength - 5;
          break;
        case 8:
          svbp.sv_pos2 = start;
          break;
        default:
          svbp.sv_pos1 = svbp.sv_pos2;
          svbp.sv_pos2 = start;
      }

      svbp.sv_chr1 = svbp.sv_chr2;
      svbp.sv_chr2 = chr;
      String sv_ori2_tmp1 = "";
      String sv_ori2_tmp2 = "";
      if (ori_readcount.containsKey("+")) sv_ori2_tmp1 = ori_readcount.get("+") + "";
      if (ori_readcount.containsKey("-")) sv_ori2_tmp2 = ori_readcount.get("-") + "";
      svbp.sv_ori2 = sv_ori2_tmp1 + "+" + sv_ori2_tmp2 + "-";

      // add up the read number
      for (Integer i_node = svbp.firstNode; i_node < node; i_node++) {
        HashMap<String, Integer> read_count_ROI_map_second =
            an.ap.rcInfo.read_count_roi_map.get(i_node);
        for (String rg : read_count_ROI_map_second.keySet()) {
          if (!svbp.read_count.containsKey(rg)) svbp.read_count.put(rg, (Integer) 0);
          svbp.read_count.put(rg, svbp.read_count.get(rg) + read_count_ROI_map_second.get(rg));
        }
        if (i_node == svbp.firstNode) continue;
        HashMap<String, Integer> read_count_FR_map_second =
            an.ap.rcInfo.read_cound_fr_map.get(i_node);
        for (String rg : read_count_FR_map_second.keySet()) {
          if (!svbp.read_count.containsKey(rg)) svbp.read_count.put(rg, (Integer) 0);
          svbp.read_count.put(rg, svbp.read_count.get(rg) + read_count_FR_map_second.get(rg));
        }

      }


    } else {
      svbp.firstNode = node;
      svbp.sv_chr1 = chr;
      svbp.sv_chr2 = chr;
      svbp.sv_pos1 = start;
      svbp.sv_pos2 = end;
      String sv_ori1_tmp1 = "";
      String sv_ori1_tmp2 = "";
      if (ori_readcount.containsKey("+")) sv_ori1_tmp1 = ori_readcount.get("+") + "";
      if (ori_readcount.containsKey("-")) sv_ori1_tmp2 = ori_readcount.get("-") + "";
      svbp.sv_ori1 = sv_ori1_tmp1 + "+" + sv_ori1_tmp2 + ("-");
      svbp.normal_rp = (int) nrp;
    }

  }

  private int calculateMaxScore() {
    double maxscore = 0;
    int flag = 0;

    for (int SVType : SVTypeReadcount.keySet()) {
      double ptype = SVTypeReadcount.get(SVType) / (double) nread_pairs;
      if (maxscore < ptype) {
        maxscore = ptype;
        flag = SVType;
      }
    }
    return flag;

  }

  private void cleanNonSupportives(Integer node) {
    ArrayList<SAMRecord> nonsupportives = new ArrayList<SAMRecord>();
    for (SAMRecord record : an.ap.regionInfo.regs.get(node))
      if (readPair.containsKey(record.getReadName())) nonsupportives.add(record);
    an.ap.regionInfo.regs.put(node, nonsupportives);
  }

  private void examinRegion(Integer node) {
    HashMap<String, Integer> orient_count = new HashMap<String, Integer>();
    ArrayList<SAMRecord> nonsupportives = new ArrayList<SAMRecord>();
    if (an.ap.regionInfo.regs.containsKey(node))
      for (SAMRecord record : an.ap.regionInfo.regs.get(node)) {
        updateNonSupportingReads(record, orient_count, nonsupportives);
      }
    an.ap.regionInfo.regs.put(node, nonsupportives);
    SVTypeOrientCounts.add(orient_count);


  }

  private void updateNonSupportingReads(SAMRecord record, HashMap<String, Integer> orient_count,
      ArrayList<SAMRecord> nonsupportives) {
    String rg = record.getReadGroup().getId();
    String ori = "+";
    if (record.getReadNegativeStrandFlag()) ori = "-";
    if (!orient_count.containsKey(ori)) orient_count.put(ori, 0);
    orient_count.put(ori, orient_count.get(ori) + 1);
    if (!readPair.containsKey(record.getReadName())) {
      readPair.put(record.getReadName(), record);
      nonsupportives.add(record);
    } else {
      int flag = (Integer) (record.getAttribute("fz"));
      if (!SVTypeReadcount.containsKey(flag)) SVTypeReadcount.put(flag, 0);
      SVTypeReadcount.put(flag, SVTypeReadcount.get(flag) + 1);

      if (!SVTypeLibraryReadcount.containsKey(flag))
        SVTypeLibraryReadcount.put(flag, new HashMap<String, Integer>());
      HashMap<String, Integer> SVTypeLibraryReadcountFlag = SVTypeLibraryReadcount.get(flag);
      if (!SVTypeLibraryReadcountFlag.containsKey(rg)) SVTypeLibraryReadcountFlag.put(rg, 0);
      SVTypeLibraryReadcountFlag.put(rg, SVTypeLibraryReadcountFlag.get(rg) + 1);

      if (!SVTypeLibraryMeanspan.containsKey(flag))
        SVTypeLibraryMeanspan.put(flag, new HashMap<String, Integer>());
      HashMap<String, Integer> SVTypeLibraryMeanspanFlag = SVTypeLibraryMeanspan.get(flag);
      if (!SVTypeLibraryMeanspanFlag.containsKey(rg)) SVTypeLibraryMeanspanFlag.put(rg, 0);
      SVTypeLibraryMeanspanFlag.put(rg,
          SVTypeLibraryMeanspanFlag.get(rg) + Math.abs(record.getInferredInsertSize()));

      nread_pairs++;
      freeReads.add(record.getReadName());
      supportReads.add(record);
      supportReads.add(readPair.get(record.getReadName()));
      readPair.remove(record.getReadName());

    }

  }
}
