/*
 * The BSD License
 * 
 * Copyright (C) 2012 Bina Technologies Inc.
 * 
 * www.binatechnologies.com
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all copies or
 * substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
 * KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

package bina.seqalto.utilities.bamextractor;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMRecordUtil;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * Reads the content of a BAM file and emits reads that have a BWA style trimmed size of > min and <
 * max.
 * 
 * @author amir@binatechnologies.com (Amirhossein Kiani)
 */
public class BamToFastQ {
  // define command line options and arguments

  @Option(name = "--in_bam_file", usage = "BAM file to extract reads from.", metaVar = "<bam_file>", required = true)
  File inBamFile;

  @Option(name = "--out_file_prefix", usage = "Output prefix.", required = true)
  String outPrefix;

  /** main function */
  public static void main(String[] args) throws IOException {
    new BamToFastQ().doMain(args);
  }

  private HashMap<String, SAMRecord> memoryPool = new HashMap<String, SAMRecord>();


  /** processes the passed command line options and runs a mode */
  private void doMain(String[] args) throws IOException {
    CmdLineParser parser = new CmdLineParser(this);
    parser.setUsageWidth(80);
    try {
      parser.parseArgument(args);

      System.out.println("In file: " + inBamFile);

      File outFile1 = new File(outPrefix + "_1.fq");
      File outFile2 = new File(outPrefix + "_2.fq");
      File discordantFile = new File(outPrefix + "_discordant_or_unmapped.bam");
      BAMFileWriter discordant = new BAMFileWriter(discordantFile);

      System.out.println("Out files: " + outFile1 + ", " + outFile2);

      checkFileExists(inBamFile);

      SAMFileReader reader = new SAMFileReader(inBamFile);

      SAMFileHeader header = reader.getFileHeader();
      discordant.setSortOrder(SortOrder.queryname, false);
      header.setSortOrder(SortOrder.queryname);
      discordant.setHeader(header);


      int pairsOutputted = 0;
      int discordantOrUnmapped = 0;

      int totalReads = 0;

      BufferedWriter writer1 = new BufferedWriter(new FileWriter(outFile1));
      BufferedWriter writer2 = new BufferedWriter(new FileWriter(outFile2));

      reader.setValidationStringency(ValidationStringency.LENIENT);



      for (SAMRecord r : reader) {
        totalReads++;
        if (r.getReadUnmappedFlag() || r.getMateUnmappedFlag()) {
          // if it's discordant or unmapped
          discordant.addAlignment(r);
          discordantOrUnmapped++;
        } else {
          if (memoryPool.containsKey(r.getReadName())) {
            // have already seen both pairs
            if (r.getSecondOfPairFlag()) {
              // if second pair is seen after first pair
              writer1.write(readToFastQ(memoryPool.get(r.getReadName())));
              writer2.write(readToFastQ(r));
            } else {
              // if second pair is seen before first pair
              writer1.write(readToFastQ(r));
              writer2.write(readToFastQ(memoryPool.get(r.getReadName())));
            }
            pairsOutputted++;
            memoryPool.remove(r.getReadName());
          } else {
            // haven't seen the read yet
            memoryPool.put(r.getReadName(), r);
          }
        }
        
        if(totalReads % 1000000 == 0){
          System.out.println("Size of memory pool is " + memoryPool.size() + " reads...");
        }
        
      }


      discordant.close();

      System.out.println("Total orphan reads (not outputted): "  + memoryPool.size() );
      System.out.println("Pairs outputted in first pass: " + pairsOutputted);
      System.out.println("Discordant or unmapped reads: " + discordantOrUnmapped);
      System.out.println("Sorted by Query name discordant reads written to: " + discordantFile);
      System.out.println("Starting to read the sorted file and appending to files!");
      SAMFileReader sortedReader = new SAMFileReader(discordantFile);

      SAMRecordIterator iter = sortedReader.iterator();
      SAMRecord firstEnd;
      SAMRecord secondEnd;
      int spacedPairs = 0;
      
      while (iter.hasNext()) {
        firstEnd = iter.next();
        try {
          secondEnd = iter.next();
        } catch (Exception e) {
          System.err.println("Orphan read found in input: " + firstEnd);
          System.err.println("Need to break :( But FASTQ data should be in good shape so far.");
          writer1.close();
          writer2.close();
          break;
        }
        if (firstEnd.getReadName().equals(secondEnd.getReadName())) {
          writer1.write(readToFastQ(firstEnd));
          writer2.write(readToFastQ(secondEnd));
          spacedPairs ++;
        } else {
          System.err.println("Orphan read found in input: " + firstEnd);
          System.err.println("Need to break :( But FASTQ data should be in good shape so far.");
          writer1.close();
          writer2.close();
          break;
        }
      }

      System.out.println("Spaced pairs were appended to FASTQ in second pass: " + spacedPairs);
      System.out.println();
      System.out.println("Total reads read from BAM: " + totalReads);
      System.out.println("Total reads outputted to FASTQ: " + (spacedPairs * 2 + pairsOutputted * 2));
      writer1.close();
      writer2.close();
    } catch (CmdLineException e) {
      parser.printUsage(System.err);
      return;
    }
  }


  private String readToFastQ(SAMRecord r) {
    if (r.getReadNegativeStrandFlag()) {
      SAMRecordUtil.reverseComplement(r);
    }
    return "@" + r.getReadName() + "/" + (r.getFirstOfPairFlag() ? "1" : "2") + "\n"
        + r.getReadString() + "\n+\n" + r.getBaseQualityString() + "\n";
  }

  private void checkFileExists(File file) {
    if (!file.exists()) {
      throw new RuntimeException("File " + file + " does not exist!");
    }
  }

}
