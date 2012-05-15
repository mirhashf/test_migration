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
import java.util.ArrayList;

import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * Reads the content of a BAM file and emits reads that have a BWA style trimmed size of > min and <
 * max.
 * 
 * @author amir@binatechnologies.com (Amirhossein Kiani)
 */
public class BamExtractor {
  // define command line options and arguments

  @Option(name = "--in_bam_file", usage = "BAM file to extract reads from.", metaVar = "<bam_file>", required = true)
  File inBamFile;

  @Option(name = "--out_bam_file", usage = "BAM file to write reads to.", metaVar = "<bam_file>", required = true)
  File outBamFile;

  @Option(name = "--max_length", usage = "Reads with trimmed-lenght of <= this length will be written to the output file. Default: 150", required = false)
  int maxLenght = 150;

  @Option(name = "--min_length", usage = "Reads with trimmed-lenght of >= this length will be written to the output file. Defualt: 0", required = false)
  int minLenght = 0;
  
  @Option(name = "--trim-quality", usage = "Quality with which to trim reads. Default: 30", required = false)
  int trimQuality = 30;

  @Option(name = "--remove-flags", usage = "Remove RG, PL, PU, LB, SM flags. Default: false", required = false)
  boolean removeFlags = false;

  
  public static String newline = System.getProperty("line.separator");

  /** main function */
  public static void main(String[] args) throws IOException {
    new BamExtractor().doMain(args);
  }

  /** processes the passed command line options and runs a mode */
  private void doMain(String[] args) throws IOException {
    CmdLineParser parser = new CmdLineParser(this);
    parser.setUsageWidth(80);
    try {
      parser.parseArgument(args);
      
      System.out.println("In file: " + inBamFile);
      System.out.println("Out file: " + outBamFile);
      System.out.println("Min Read Length: " + minLenght);
      System.out.println("Max Read Length: " + maxLenght);
      System.out.println("Trim quality: " + trimQuality);
      
      checkFileExists(inBamFile);
      
      SAMFileReader reader = new SAMFileReader(inBamFile);
      BAMFileWriter writer = new BAMFileWriter(outBamFile);
      
      reader.setValidationStringency(ValidationStringency.LENIENT);
      
      SAMFileHeader header = reader.getFileHeader().clone();
      
      if(removeFlags){
       header.setReadGroups(new ArrayList<SAMReadGroupRecord>());
      }
      
      writer.setHeader(header);
      
      for(SAMRecord r : reader){
        int trimmedLength = getBwaTrimmedSize(trimQuality, r);
        if(trimmedLength >= minLenght && trimmedLength <= maxLenght){
          if(removeFlags){
            r.setAttribute("RG", null);
            r.setAttribute("PL", null);
            r.setAttribute("PU", null);
            r.setAttribute("SM", null);
          }
          writer.addAlignment(r);
        }
      }
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
