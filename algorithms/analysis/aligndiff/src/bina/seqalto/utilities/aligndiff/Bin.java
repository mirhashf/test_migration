/*
 * The BSD License
 * 
 * Copyright (C) 2012 Bina Technologies Inc. www.binatechnologies.com
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

package bina.seqalto.utilities.aligndiff;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

/**
 * Accumulates results for each "Bin" of SAM records and calculates statistics
 * on mapping quality, number of unmapped, etc.
 * 
 * @author amir@binatechnologies.com (Amirhossein Kiani)
 */
public class Bin {
  // buffer size for flushing
  private final int BUFFER_SIZE = 4000;
  private final String nl = System.getProperty("line.separator");
  public long CountReads = 0;
  public int discordantReads = 0;
  public long SumMappingQuality = 0;
  public long SumBaseQuality = 0;
  public long CountBases = 0;
  public long CountUnmapped = 0;
  public long CountMultiplyMapped = 0;
  public long SumGapOpen = 0;
  public long SumEditDistance = 0;
  public long SumMismatches = 0;
  public long SumMostRepetetiveKmers = 0;
  private boolean isPersistent = false;
  private BufferedWriter fileWriter = null;
  private String filePath = null;
  public long SumGapExtensions = 0;
  public Bin lowMapqBin;
  private SAMFileReader samReader = null;
  public Bin shortReadsBin;
  public Bin unmappedReadsBin;
  private boolean isMainBin = false;
  public long SumSoftClipping = 0;
  public double total = 0;
  private boolean outputPair = false;

  public double SumNormalizedEditDistance = 0;

  /* constructor for non-persistent bin (a bin that doesn't write output to a file) */
  public Bin() {

  }

  /* constructor for persistent/non-persistent bin */
  public Bin(String binFilePath, boolean mainBin, SAMFileReader reader, boolean outputPair)
      throws IOException {
    
    this.isPersistent = true;
    this.fileWriter = new BufferedWriter(new FileWriter(binFilePath));
    this.filePath = binFilePath;
    this.outputPair = outputPair;
    if (mainBin) {
      isMainBin = true;
      lowMapqBin = new Bin(binFilePath + "_LowMappingQuality", false, reader, outputPair);
      shortReadsBin = new Bin(binFilePath + "_ShortReads", false, reader, outputPair);
      unmappedReadsBin = new Bin(binFilePath + "_UnmappedReads", false, reader, outputPair);
    }
    samReader = reader;
  }

  public void close() throws IOException {
    if (this.fileWriter != null) {
      this.fileWriter.flush();
      this.fileWriter.close();
      this.lowMapqBin.fileWriter.flush();
      this.lowMapqBin.fileWriter.close();
    }
  }

  /* Inserts a SAM record to this @Bin */ 
  public void insert(SAMRecord record) throws IOException {
    total++;
    if (Math.abs(record.getAlignmentStart() - record.getMateAlignmentStart()) > 500) {
      discordantReads++;
    }

    if (record.getReadUnmappedFlag() && isMainBin) {
      this.unmappedReadsBin.insert(record);
      return;
    }

    if (record.getReadLength() < -1 && isMainBin) {
      this.shortReadsBin.insert(record);
      return;
    }

    if (record.getMappingQuality() < 5 && isMainBin) {
      this.lowMapqBin.insert(record);
      return;
    }

    if (isPersistent) {
      this.fileWriter.write(record.getSAMString());
      if (outputPair) {
        SAMRecord mate = samReader.queryMate(record);
        if (mate != null) {
          this.fileWriter.write(mate.getSAMString());
        } else {
          System.err.println("Could not find mate for " + record.getReadName() + " in " + filePath);
        }
      }
      if (CountReads % BUFFER_SIZE == 0) {
        fileWriter.flush();
      }
    }

    CountReads++;
    if (!record.getReadUnmappedFlag()) {
      // Sam flag indicates mapped
      SumGapOpen += SAMUtil.getIntegerFlagSafe(record, "XO");
      SumEditDistance += SAMUtil.getIntegerFlagSafe(record, "NM");
      SumMismatches += SAMUtil.getIntegerFlagSafe(record, "XM");
      SumGapExtensions += SAMUtil.getIntegerFlagSafe(record, "XG");
      SumMappingQuality += record.getMappingQuality();
      double clippingSize = getClippingSize(record);
      SumSoftClipping += clippingSize;
      SumNormalizedEditDistance +=
          (SAMUtil.getIntegerFlagSafe(record, "NM") / (record.getReadLength() - clippingSize));
    }
    SumMostRepetetiveKmers += calculateRepetetiveness(record.getReadString());

    // both for mapped and unmapped reads

    // assuming Sanger base qualities
    char[] baseQualities = record.getBaseQualityString().toCharArray();

    for (char baseQuality : baseQualities) {

      if (baseQuality != '#') {
        SumBaseQuality += (baseQuality - 33);
        CountBases++;
      }
    }

  }

  /* Calculates repetitiveness measure for a string of bases 
   * This is done by breaking the bases into overlapping kmers of size 5 and then
   * counting the occurrences of each kmer.
   */
  public static int calculateRepetetiveness(String bases) {
    HashMap<String, Integer> kmerCountMap = new HashMap<String, Integer>();
    for (int i = 0; i <= bases.length() - 5; i++) {
      String kmer = bases.substring(i, i + 5);
      if (kmerCountMap.containsKey(kmer)) {
        kmerCountMap.put(kmer, kmerCountMap.get(kmer) + 1);
      } else {
        kmerCountMap.put(kmer, 1);
      }
    }
    int max = 0;
    // find max
    for (Integer val : kmerCountMap.values()) {
      if (val > 2) {
        max += val;
      }
    }
    return max;
  }

  public String safeAverage(double nominator, double denominator) {
    DecimalFormat formatter = new DecimalFormat("###.###");
    if (denominator != 0) {
      return String.valueOf(formatter.format(nominator / denominator));
    } else {
      return "Division by zero! (" + nominator + "/0)";
    }
  }

  private double getClippingSize(SAMRecord rec) {
    double unclippedStart = rec.getUnclippedStart();
    double unclippedEnd = rec.getUnclippedEnd();
    double start = rec.getAlignmentStart();
    double end = rec.getAlignmentEnd();
    if (unclippedStart == start && unclippedEnd == end) {
      return 0;
    } else {
      if (start != unclippedStart) {
        return Math.abs(unclippedStart - start);
      } else {
        return Math.abs(end - unclippedEnd);
      }
    }
  }

  /* Returns statistics for this bin as a String*/
  public String toString(double total) {
    DecimalFormat formatter = new DecimalFormat("###.###");

    String out = "";
    if (isMainBin) {
      out += "==================================";

      if (filePath != null) {
        out += filePath;
      }
      out += "==================================" + nl;
    }
    out += "Bin Size:" + this.total / total * 100 + "%" + nl;
    out +=
        "CountReads: " + this.CountReads + "(" + formatter.format(CountReads / total * 100.) + "%)"
            + nl;
    out +=
        "CountUnmapped: " + this.CountUnmapped + "("
            + formatter.format(CountUnmapped / total * 100.) + "%)" + nl;
    out +=
        "Sum Gap Opens (xo): " + SumGapOpen + "(Average over mapped reads = "
            + safeAverage(SumGapOpen, (CountReads - CountUnmapped)) + ")" + nl;
    out +=
        "Sum Edit Distance (nm): " + SumEditDistance + "(Average over mapped reads = "
            + safeAverage(SumEditDistance, (CountReads - CountUnmapped)) + ")" + nl;
    out +=
        "Sum Mistmach (xm): " + SumMismatches + "(Average over mapped reads = "
            + safeAverage(SumMismatches, (CountReads - CountUnmapped)) + ")" + nl;

    out +=
        "Gap Extensions (xg): " + SumGapExtensions + "(Average over mapped reads = "
            + safeAverage(SumGapExtensions, (CountReads - CountUnmapped)) + ")" + nl;
    out +=
        "Sum Base Qualities: " + SumBaseQuality + "(Average over all reads = "
            + safeAverage(SumBaseQuality, CountBases) + ")" + nl;

    out +=
        "Sum Mapping Quality: " + SumMappingQuality + "(Average over all reads = "
            + safeAverage(SumMappingQuality, CountReads) + ")" + nl;

    out +=
        "Sum Soft Clipping: " + SumSoftClipping + "(Average over all reads = "
            + safeAverage(SumSoftClipping, CountReads - CountUnmapped) + ")" + nl;

    out += "Average read length: " + safeAverage(CountBases, CountReads - CountUnmapped) + nl + nl;

    out += "Discordant Reads: " + discordantReads + nl;

    if (isMainBin) {

      out += ">>>---- Low MapQ bin ----<<<" + nl;
      out += this.lowMapqBin.toString(total);
      if (this.unmappedReadsBin.total > 0) {
        out += ">>>---- Unmapped reads bin ----<<<" + nl;
        out += this.unmappedReadsBin.toString(total);
      }

    }

    return out;
  }
}
