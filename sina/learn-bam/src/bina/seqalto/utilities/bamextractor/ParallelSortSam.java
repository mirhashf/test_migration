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
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.PriorityQueue;

import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordCoordinateComparator;
import net.sf.samtools.SAMRecordQueryNameComparator;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * 
 * @author amir@binatechnologies.com (Amirhossein Kiani)
 */
public class ParallelSortSam {
  // define command line options and arguments

  @Option(name = "--in_bam_file", usage = "BAM file to sort.", metaVar = "<bam_file>", required = true)
  File inBamFile;

  @Option(name = "--out_bam_file", usage = "BAM file to write to.", metaVar = "<bam_file>", required = true)
  File outBamFile;

  @Option(name = "--threads", usage = "Number of threads (default: 8).", metaVar = "<bam_file>", required = false)
  int numThreads = 8;

  @Option(name = "--sort-order", usage = "Sort order: coordinate, queryname (default: coordinate)", required = false)
  String sortOrder = "coordinate";

  @Option(name = "--max-ram-records", usage = "Total number of records in ram.", required = false)
  int maxRecordsInRam = 1000000;


  /**
   * main function
   * 
   * @throws InterruptedException
   */
  public static void main(String[] args) throws IOException, InterruptedException {
    new ParallelSortSam().doMain(args);
  }

  /**
   * processes the passed command line options and runs a mode
   * 
   * @throws InterruptedException
   */
  private void doMain(String[] args) throws IOException, InterruptedException {
    CmdLineParser parser = new CmdLineParser(this);
    parser.setUsageWidth(80);
    try {
      parser.parseArgument(args);
      checkFileExists(inBamFile);
      ArrayList<SorterConsumer> consumers = new ArrayList<ParallelSortSam.SorterConsumer>();


      SAMFileReader reader = new SAMFileReader(inBamFile);
      SAMFileHeader header = reader.getFileHeader();
      Comparator<SAMRecord> comparator;

      if (sortOrder.equals("coordinate")) {
        comparator = new SAMRecordCoordinateComparator();
      } else if (sortOrder.equals("queryname")) {
        comparator = new SAMRecordQueryNameComparator();
      } else {
        System.err.println("Invalid sort order!");
        return;
      }


      Calendar now = Calendar.getInstance();


      System.out.println(now.getTime() + " Started reading input.");

      // start threads that read
      for (int i = 0; i < numThreads; i++) {
        SorterConsumer consumer = new SorterConsumer(comparator);
        consumers.add(consumer);
      }

      int packet = 0;

      for (SAMRecord read : reader) {
        consumers.get(packet % numThreads).sortedList.add(read);
      }

      for (SorterConsumer consumer : consumers) {
        consumer.start();
      }

      for (SorterConsumer consumer : consumers) {
        consumer.join();
      }

      

      System.out.println(now.getTime() + " Finished reading the file in...");

      BAMFileWriter merger = new BAMFileWriter(outBamFile);

      merger.setHeader(header);
      merger.setSortOrder(header.getSortOrder(), true);

      PriorityQueue<SAMRecord> mergedHeap =
          new PriorityQueue<SAMRecord>(maxRecordsInRam, comparator);

      for (SorterConsumer consumer : consumers) {
        while (!consumer.sortedList.isEmpty()) {
          mergedHeap.add(consumer.sortedList.remove(0));
        }
      }

      System.out.println(now.getTime() + " Finished merging...");

      // write out
      for (SAMRecord read : mergedHeap) {
        merger.addAlignment(read);
      }

      System.out.println(now.getTime() + " Finished writign...");


      merger.close();
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


  class SorterConsumer extends Thread {
    private ArrayList<SAMRecord> sortedList;
    private Comparator<SAMRecord> comparator;

    SorterConsumer(Comparator<SAMRecord> comparator) {
      this.comparator = comparator;
      this.sortedList = new ArrayList<SAMRecord>(maxRecordsInRam / maxRecordsInRam);
    }

    @Override
    public void run() {
      Collections.sort(sortedList, comparator);
    }

  }
}
