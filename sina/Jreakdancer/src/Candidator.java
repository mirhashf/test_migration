import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.security.AllPermission;
import java.util.*;

import javax.swing.SortOrder;

import net.sf.picard.sam.SamPairUtil;
import net.sf.picard.sam.SamPairUtil.PairOrientation;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.BAMIndexer;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;

public class Candidator {

  private static final boolean VERBOSE=true;
  
  public static final Integer numReadsToProcess = 20000;
 
  
  private SAMFileHeader header;
 
  private CandidatorParameters cp=new CandidatorParameters();

  protected static final String CANDIDATE_PARAMETER_FILE="candidate_parameters.obj";
  
  public Candidator(String inputBAM, String inputStats, String scratchDir, String chr, Integer start,
      Integer end, int cut_sd, int minMapQ) {
    super();
    cp.cut_sd = cut_sd;
    cp.MINMAPQ = minMapQ;
    cp.inputStats = inputStats;
    cp.inputBAM = inputBAM;
    cp.chr = chr;
    cp.start = start;
    cp.end = end;
    cp.scratchDir = scratchDir;
    
    getChrLengths();
    cp.refLen = cp.chrLength.get(chr);
   
    File scratchDirFile = new File(scratchDir);
    if (!scratchDirFile.exists()) scratchDirFile.mkdirs();
    
    cp.readgroup = new HashMap<String, ReadGroupInfo>();
    extractReadGroupStatistics();
    extractCandidateReads();
  }

  public static void getBamStatistics(String[] args) throws IOException {
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
    Candidator candidator =
        new Candidator(inputBam, inputStats, scratchDir, chr, start, end, cut_sd, minMapQ);

  }


  private void extractReadGroupStatistics() {
    try {
      BufferedReader r = new BufferedReader(new FileReader(cp.inputStats));
      String line = "";
      while ((line = r.readLine()) != null) {
        String[] tokens = line.split("\\t");
        if (tokens.length < 4) continue;
        String rg = tokens[0];
        Double mean = Double.valueOf(tokens[1]);
        Double sd = Double.valueOf(tokens[2]);
        Double readLen = Double.valueOf(tokens[3]);
        ReadGroupInfo info = new ReadGroupInfo(0);
        info.mean = mean;
        info.std = sd;
        info.readlen = readLen;
        info.upper = info.mean + info.std * cp.cut_sd;
        info.lower = info.mean - info.std * cp.cut_sd;
        if (info.lower < 0) info.lower = 0;
        cp.readgroup.put(rg, info);
        if (cp.chrDistanceToSeeAnomaly > info.mean - 2 * info.readlen)
          cp.chrDistanceToSeeAnomaly = info.mean - 2 * info.readlen;

      }
      r.close();
      if (cp.chrDistanceToSeeAnomaly < 50) cp.chrDistanceToSeeAnomaly = 50;
    } catch (IOException e) {
      e.printStackTrace();
      System.exit(1);
    }

  }

  private void extractCandidateReads() {
     
    SAMFileReader reader = new SAMFileReader(new File(cp.inputBAM));
    SAMRecordIterator iterator = reader.iterator();
    if (cp.chr.trim().length() > 0) {
      if (!reader.hasIndex())
        BAMIndexer.createAndWriteIndex(new File(cp.inputBAM), new File(cp.inputBAM + "bai"), false);
      iterator.close();
      iterator = reader.query(cp.chr, (int) cp.start, (int) cp.end, true);
    }

    while (iterator.hasNext()) {
      SAMRecord record = iterator.next();
      int flag=AlnParser.alnParse(record);
      String rg=record.getReadGroup().getId();
      
      if (record.getMappingQuality() >= cp.MINMAPQ && flag < 32 && flag >= 18) {
        if(!cp.numAllReadsPerLibrary.containsKey(rg))
          cp.numAllReadsPerLibrary.put(rg, (Integer)0);
       cp.numAllReadsPerLibrary.put(rg, cp.numAllReadsPerLibrary.get(rg)+1);
       continue;
      }
      if(record.getMappingQuality()<cp.MINMAPQ || record.getDuplicateReadFlag()==true || flag==0)
        continue;
  
      if(record.getReadUnmappedFlag() || record.getMateUnmappedFlag() || !record.getReferenceName().equals(record.getReferenceName()))
        continue;

      int insertSize = Math.abs(record.getInferredInsertSize());
      if(!cp.readgroup.containsKey(rg)) continue;
      
      ReadGroupInfo info = cp.readgroup.get(rg);
      if(insertSize>info.upper && flag==18)
        flag=2;
      if(insertSize<= info.upper && flag==2)
        flag=18;
      
      if(insertSize<info.lower && flag==18)
        flag=3;
      if(insertSize>= info.lower && flag==3)
        flag=18;
      if(flag==20)
        flag=4;
  
      if(flag==18)
        continue;

      if(!cp.numAllBadReadsPerType.containsKey(flag))
        cp.numAllBadReadsPerType.put(flag, new HashMap<String,Integer>());
      if(!cp.numAllBadReadsPerType.get(flag).containsKey(rg))
        cp.numAllBadReadsPerType.get(flag).put(rg, (Integer)0);
      cp.numAllBadReadsPerType.get(flag).put(rg, cp.numAllBadReadsPerType.get(flag).get(rg)+1);
      
      cp.numReads++;
      if (cp.numReads > numReadsToProcess) {
        cp.refLen = record.getAlignmentStart() - cp.start;
        System.out.println("refLen:" + cp.refLen);
        break;
      }
    }
    
    iterator.close();
    reader.close();

    for(String rg:cp.numAllReadsPerLibrary.keySet()){
      cp.readDensity.put(rg, (double)cp.numAllReadsPerLibrary.get(rg)/cp.refLen);
      if(VERBOSE==true){
        double seq_coverage=cp.readDensity.get(rg)*cp.readgroup.get(rg).readlen;
        double physical_coverage=cp.readDensity.get(rg)*cp.readgroup.get(rg).mean;
        System.err.println(rg+"\t"+seq_coverage+"\t"+physical_coverage);
        for(Integer svType:cp.numAllBadReadsPerType.keySet())
          System.err.println(svType+"\t"+rg+"\t"+cp.numAllBadReadsPerType.get(svType).get(rg));
      }
     
      int numReadsDiscrepantSize=0;
      try{
        numReadsDiscrepantSize+=cp.numAllBadReadsPerType.get(2).get(rg)+cp.numAllBadReadsPerType.get(3).get(rg);
        double discrepantDensity=(double)cp.refLen/numReadsDiscrepantSize;
        if(cp.chrDistanceToSeeAnomaly>discrepantDensity )
          cp.chrDistanceToSeeAnomaly=discrepantDensity;
      }catch(Exception e){}
    }
    
    
    try {
      writeCandidates();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  private void writeCandidates() throws IOException {
   
    ObjectOutputStream w2 =
       new ObjectOutputStream(new FileOutputStream(
            new File(cp.scratchDir, CANDIDATE_PARAMETER_FILE).getAbsolutePath()));
    w2.writeObject(cp);
    w2.close();
  }

  private void getChrLengths() {
    cp.chrLength = new HashMap<String, Integer>();
    SAMFileReader reader = new SAMFileReader(new File(cp.inputBAM));
    header = reader.getFileHeader();

    if (cp.chr.equals("")) {
      for (SAMSequenceRecord chr : header.getSequenceDictionary().getSequences())
        cp.chrLength.put(chr.getSequenceName(), (Integer) chr.getSequenceLength());

    } else {
      if (cp.start != 0 || cp.end != 0)
        cp.chrLength.put(cp.chr, cp.end - cp.start);
      else
        cp.chrLength.put(cp.chr, (Integer) header.getSequence(cp.chr).getSequenceLength());

    }
    reader.close();
    System.out.println(cp.chrLength);
  }


}


