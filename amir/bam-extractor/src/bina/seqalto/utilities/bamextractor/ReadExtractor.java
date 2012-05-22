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

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * Reads the content of a BAM file and emits reads that have a BWA style trimmed size of > min and <
 * max.
 * 
 * @author amir@binatechnologies.com (Amirhossein Kiani)
 */
public class ReadExtractor {
  // define command line options and arguments

  @Option(name = "--in_bam_file", usage = "BAM file to extract reads from.", metaVar = "<bam_file>", required = true)
  File inBamFile;

  @Option(name = "--start", usage = "Start pos", required = true)
  int start;

  @Option(name = "--end", usage = "End pos", required = true)
  int end;

  @Option(name = "--ref", usage = "Reference", required = true)
  String ref;


  public static String newline = System.getProperty("line.separator");

  /** main function */
  public static void main(String[] args) throws IOException {
    new ReadExtractor().doMain(args);
  }

  /** processes the passed command line options and runs a mode */
  private void doMain(String[] args) throws IOException {
    CmdLineParser parser = new CmdLineParser(this);
    parser.setUsageWidth(80);
    try {
      parser.parseArgument(args);
      checkFileExists(inBamFile);
      SAMFileReader reader = new SAMFileReader(inBamFile);
      reader.setValidationStringency(ValidationStringency.LENIENT);

      SAMRecordIterator results = reader.queryOverlapping(ref, start, end);
      SAMRecord r;
      while (results.hasNext()) {
        r = results.next();
        System.out.println(r.getReadString());

      }


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
}
