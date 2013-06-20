import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import net.sf.samtools.BAMIndex;
import net.sf.samtools.BAMIndexer;
import net.sf.samtools.BrowseableBAMIndex;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import java.util.*;

public class BamStatistics {

  private int MINMAPQ = 50;
  private static final int numReadsToProcess = 1000000;

  private HashMap<String, ReadGroupInfo> read_group;
  private static Random rand = new Random();

  protected String chr;
  protected Integer start;
  protected Integer end;
  protected String inputBam;
  protected PrintWriter outputFile;


  public BamStatistics(String chr, Integer start, Integer end, String inputBam, PrintWriter outputFile) {
    super();
    this.chr = chr;
    this.start = start;
    this.end = end;
    this.inputBam = inputBam;
    this.outputFile = outputFile;

    read_group = new HashMap<String, ReadGroupInfo>();

    SAMFileReader reader = new SAMFileReader(new File(inputBam));
    SAMRecordIterator iterator = reader.iterator();
    if (chr.trim().length() > 0) {
      if (!reader.hasIndex())
        BAMIndexer.createAndWriteIndex(new File(inputBam), new File(inputBam + "bai"), false);
      iterator.close();
      iterator = reader.query(chr, (int) start, (int) end, true);
    }
    int count = 0;
    while (iterator.hasNext()) {
      if (rand.nextDouble() > 0.1) continue;
      count++;
      if (count > numReadsToProcess) break;
      SAMRecord record = iterator.next();
      String rg = record.getReadGroup().getReadGroupId();
      if (!read_group.containsKey(rg)) read_group.put(rg, new ReadGroupInfo(MINMAPQ));
      ReadGroupInfo info = read_group.get(rg);
      info.addInfo(record);
    }

    for (String rg : read_group.keySet()) {
      ReadGroupInfo info = read_group.get(rg);
      info.mean /= info.count;
      info.std = Math.sqrt(info.std / info.count - info.mean * info.mean);
      info.readlen /= info.count;
    }

    for (String rg : read_group.keySet()) {
      ReadGroupInfo info = read_group.get(rg);
      writeOutput(rg,info,outputFile);
      System.out.println(rg + "\t" + info);
    }
    
    reader.close();
   
  }


  private void writeOutput(String rg, ReadGroupInfo info, PrintWriter outputFile) {
    System.err.println(rg+"\t"+info.mean+"\t"+info.std+"\t"+info.readlen);
    outputFile.println(rg+"\t"+info.mean+"\t"+info.std+"\t"+info.readlen);
  }

public static void getBamStatistics(String[] args) throws IOException {
   int count=0;
    String inputBam = args[count++];
    String output = args[count++];
    PrintWriter outputFile = new PrintWriter(System.out);
    if (!output.equals("stdout")) outputFile = new PrintWriter(new FileWriter(output));
    String scratchDir=args[count++];
    String chr = "";
    if (args.length > count) chr = args[count++];
    Integer start = 0;
    if (args.length > count) start = Integer.valueOf(args[count++]);
    Integer end = 0;
    if (args.length > count) start = Integer.valueOf(args[count]);
    BamStatistics stats = new BamStatistics(chr, start, end, inputBam, outputFile);
    outputFile.close();
    System.err.println("Bam statsitcs are generated and written to disk");
  }
}



