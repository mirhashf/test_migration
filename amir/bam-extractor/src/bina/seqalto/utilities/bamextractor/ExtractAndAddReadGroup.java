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

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * Reads the content of a BAM file and emits reads adding readgroup based on Illumina read ID
 * 
 * @author amir@binatechnologies.com (Amirhossein Kiani)
 */
public class ExtractAndAddReadGroup {
  // define command line options and arguments

  @Option(name = "--in_bam_file", usage = "BAM file to extract reads from.", metaVar = "<bam_file>", required = true)
  File inBamFile;

  @Option(name = "--out_bam_file", usage = "BAM file to write reads to.", metaVar = "<bam_file>", required = true)
  File outBamFile;

  @Option(name = "--reference", usage = "Export only reads from this reference", required = false)
  String reference;

  @Option(name = "--sample", usage = "Sample id", required = true)
  String sample;

  
  @Option(name = "--platform-unit", usage = "platform name", required = true)
  String platformUnit;


  @Option(name = "--library", usage = "library name", required = true)
  String library;

  public static String newline = System.getProperty("line.separator");

  /** main function */
  public static void main(String[] args) throws IOException {
    new ExtractAndAddReadGroup().doMain(args);
  }

  /** processes the passed command line options and runs a mode */
  private void doMain(String[] args) throws IOException {
    CmdLineParser parser = new CmdLineParser(this);
    parser.setUsageWidth(80);
    try {
      parser.parseArgument(args);
      
      Set<String> readGroups = new HashSet<String>();
      
      System.out.println("In file: " + inBamFile);
      System.out.println("Out file: " + outBamFile);
      System.out.println("Target Reference " + reference);
      
      checkFileExists(inBamFile);
      
      SAMFileReader reader = new SAMFileReader(inBamFile);
      BAMFileWriter writer = new BAMFileWriter(outBamFile);
      
      reader.setValidationStringency(ValidationStringency.LENIENT);
      
      SAMFileHeader header = reader.getFileHeader().clone();
      

      if(header.getSortOrder() == null){
        header.setSortOrder(SortOrder.unsorted);
      }

      for(int i = 1 ;i < 9; i++){
        SAMReadGroupRecord readGroupObject = new SAMReadGroupRecord(String.valueOf(i));
        readGroupObject.setSample(sample);
        readGroupObject.setLibrary(library);
        readGroupObject.setPlatformUnit(platformUnit);
        header.addReadGroup(new SAMReadGroupRecord(String.valueOf(i)));
      }

      
      writer.setHeader(header);

      for(SAMRecord r : reader){
        if((reference == null || r.getReferenceName().equals(reference)) && !r.getReadUnmappedFlag()){
          String flowcell = r.getReadName().split(":")[1];
          r.setAttribute("RG", flowcell);
          r.setAttribute("PU", platformUnit);
          r.setAttribute("SM", sample);
          r.setAttribute("LB", library);
          
          writer.addAlignment(r);
          readGroups.add(flowcell);
        }
      }
      
      System.out.println("Observed flowcells: " + readGroups);
      
      
      writer.setHeader(header);
      writer.close();
    } catch (CmdLineException e) {
      parser.printUsage(System.err);
      return;
    }
  }

  private void checkFileExists(File file) {
    if (!file.exists()) {
      throw new RuntimeException("File " + file + " does not exist!");
    }
  }

  /**
   * Get Trimed Read size then trimmed BWA Style. Inspired
   * by:http://wiki.bioinformatics.ucdavis.edu/index.php/TrimBWAstyle.pl
   */
  public static int getBwaTrimmedSize(int opt_q, SAMRecord r) {
    byte[] qs;
    int pos, maxPos, area, maxArea;
    qs = r.getBaseQualities();
    pos = qs.length;
    maxPos = pos;
    area = 0;
    maxArea = 0;
    while (pos > 0 && area >= 0) {
      area += opt_q - qs[pos-1];
      if (area > maxArea) {
        maxArea = area;
        maxPos = pos;
      }
      pos--;
    }
    return maxPos;
  }
}
