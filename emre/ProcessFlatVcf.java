import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.StringTokenizer;

public class ProcessFlatVcf {
  
  public static final String[] infoKeys = new String[]{"ABHet", "ABHom", "AC", "AF", "AN", "BaseQRankSum", "DB", "DP", 
                                                       "DS", "Dels", "FS", "HRun", "HaplotypeScore", "InbreedingCoeff", 
                                                       "MQ", "MQ0", "MQRankSum", "OND", "QD", "ReadPosRankSum", "SB"};
  
  public static void main(String[] args) throws Exception {
    String filePath = (args.length > 0 ? args[0] : null);
    if(filePath == null || filePath.trim().length() == 0) {
      throw new IllegalArgumentException("Invalid input path");
    }
    String outFilePath = (args.length > 1 ? args[1] : filePath + ".out");
    
    System.out.println("input file: " + filePath);
    System.out.println("output file: " + outFilePath);
    
    StringTokenizer lt = null;
    StringTokenizer gt = null;
    
    BufferedReader reader = null;
    BufferedWriter writer = null;
    try {
      long st = System.currentTimeMillis();
      reader = new BufferedReader(new FileReader(filePath)); 
      writer = new BufferedWriter(new FileWriter(outFilePath)); 
      
      String line;
      int i = 0;
      int headerCount = 0;
      while((line = reader.readLine()) != null) {
        if(line.startsWith("#")) {
          headerCount++;
          continue;
        }
        
        if(i > 0) {
          writer.newLine();
        }
        
        StringBuilder lineSb = null;
        lt = new StringTokenizer(line, "\t");
        while(lt.hasMoreTokens()) {
          lineSb = new StringBuilder();
          
          String chr = lt.nextToken();         
          String pos = lt.nextToken();
          String id = lt.nextToken();
          String ref = lt.nextToken();
          String alt = lt.nextToken();
          String qual = lt.nextToken();        
          String filter = lt.nextToken();
          String info = lt.nextToken();
          String format = lt.nextToken();
          String sample = lt.nextToken();
          
          gt = new StringTokenizer(info, ";");
          Map<String,String> infoMap = new HashMap<String,String>();
          while(gt.hasMoreTokens()) {
            String infoToken = gt.nextToken();
            String[] kv = infoToken.split("=");
            if(kv.length == 2) {
              infoMap.put(kv[0], kv[1]);
            } else if(kv.length == 1) {
              infoMap.put(kv[0], "1");
            }
          }
          
          StringBuilder infoSb = new StringBuilder();
          for(String infoKey : infoKeys) {
            String infoVal = infoMap.get(infoKey);
            if(infoVal == null) {
              if("DB".equals(infoKey) || "DS".equals(infoKey)) {
                infoSb.append("0");
              } else {
                infoSb.append(" "); 
              }
            } else {
              infoSb.append(infoVal);
            } 
            infoSb.append("\t");
          }
          
          gt = new StringTokenizer(sample, ":");
          StringBuilder sampleSb = new StringBuilder();
          while(gt.hasMoreTokens()) {
            String sampleToken = gt.nextToken();
            sampleSb.append(sampleToken).append("\t");
          }
          
          lineSb.append(chr).append("\t").append(pos).append("\t").append(id).append("\t").append(ref)
              .append("\t").append(alt).append("\t").append(qual).append("\t").append(filter)
              .append("\t").append(infoSb).append(sampleSb);
          
          String vcfLine = lineSb.substring(0, lineSb.length()-1);
          writer.write(vcfLine);
          
          i++;
        }
      }
      
      System.out.println("Header lines: " + headerCount);
      long et = System.currentTimeMillis();
      System.out.println("Processing completed in " + (int)Math.ceil(((et - st)/1000)) + " seconds");
    } finally {
      if(reader != null) {
        reader.close();
      }
      
      if(writer != null) {
        writer.close();
      }
    }
  }
}
