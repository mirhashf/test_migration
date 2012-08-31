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

package bina.seqalto.utilities.aligndiff;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.RuntimeEOFException;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * Basic utility to compare two BAM files read-by-read (helpful for comparing mapper's results on
 * real data) Each BAM file needs to be sorted by Read Id before the comparison is started and both
 * files need to contain the same reads (possibly with different alignments)
 * 
 * Reads are grouped into @Bin's based on different criteria and each bin calculates basic
 * statistics on it's contents.
 * 
 * Can output the results either in text or HTML (html can be activated by passing true to -output-html option)
 * 
 * @author amir@binatechnologies.com (Amirhossein Kiani)
 */
public class AlignDiff {

  // define command line options and arguments

  @Option(name = "--sort", usage = "Sort BAM file based on query name.", metaVar = "<bam_file>")
  public boolean sort = false;

  @Option(name = "--assume-sorted", usage = "Assume sorted BAMs by query name.")
  public boolean assumeSorted = false;

  @Option(name = "--diff", usage = "Diff two BAM files on different meterics by deviding them into bin's of similar characteristics. Each bin's contents is written to disk as well.", metaVar = "<bam_file_1> <bam_file_2> [output_prefix]")
  public boolean diff = false;

  @Option(name = "--stats", usage = "Stats for one BAM file")
  public boolean stats = false;

  @Option(name = "--max-in-ram", usage = "Max records in RAM.")
  public int maxRecordsInRam = 500000;

  @Option(name = "--wiggle", usage = "Number of bases to allow as wiggle room")
  public int wiggle = 10;

  @Option(name = "--output-html", usage = "Write HTML output to standard out")
  public boolean outputHTML = false;

  // uses a temp dir on home dir if it's not defined (same as what Picard does)
  @Option(name = "--temp-dir", usage = "Temporary dir")
  public File TMP_DIR = (System.getProperty("java.io.tmpdir").endsWith(
      "/" + System.getProperty("user.name"))
      ? new File(System.getProperty("java.io.tmpdir"))
      : new File(System.getProperty("java.io.tmpdir"), System.getProperty("user.name")));


  @Option(name = "--output-pairs-in-bin", usage = "Outputs Pairs for each bin corresponding the the records that do not match (as opposed to simplified CSV output). Enablign this option will slowdown the execution significantly")
  public boolean outputPair = false;

  @Option(name = "--min-read-length", usage = "Minimum read length to consider")
  public int minReadLength = 2;

  @Option(name = "--max-read-length", usage = "Maximum read length to consider")
  public int maxReadLength = 100;

  @Option(name = "--only-both-pair-mapped", usage = "Count mapped only if both pairs are mapped")
  public boolean onlyBothPair = false;

  @Option(name = "--find-discordant", usage = "Find reads that are discordant in one file and not in the other",  metaVar = "<bam_file_1> <bam_file_2> [output_prefix]")
  public boolean findDiscordant = false;
  
  @Option(name = "--max-concordance-space", usage = "Maximum distance between start positions in each read and it's pair to be considered concordant")
  public int maxConcordance = 500;

  @Argument(required = true, usage="BAM/SAM file(s) to read from")
  public List<String> arguments = new ArrayList<String>();


  // variables to keep track of global counts

  // reads in left file with MAPQ>20
  private double q20Left = 0;
  // reads in left file with MAPQ>10
  private double q10Left = 0;
  // unique reads in left file
  private double uniqueLeft = 0;
  // mapped reads in left file
  private double mappedLeft = 0;

  // reads in right file with MAPQ>20
  private double q20Right = 0;
  // reads in right file with MAPQ>10
  private double q10Right = 0;
  // unique reads in right file
  private double uniqueRight = 0;
  // mapped reads in right file
  private double mappedRight = 0;

  public static String newline = System.getProperty("line.separator");

  /** main function */
  public static void main(String[] args) throws IOException {
    new AlignDiff().doMain(args);
  }

  /** processes the passed command line options and runs a mode */
  private void doMain(String[] args) throws IOException {
    CmdLineParser parser = new CmdLineParser(this);
    parser.setUsageWidth(80);

    File fileLeft;
    try {
      // parse the arguments.
      parser.parseArgument(args);
      fileLeft = new File(arguments.get(0));
      checkFileExists(fileLeft);

      // decide what mode to run

      if (sort) {
        // sort the passed bam file
        SAMUtil.sortByQueryName(fileLeft, TMP_DIR, maxRecordsInRam);
      }

      if (stats) {
        // only output general statistics about ONE input file (#mapped, #unique, etc...)
        SAMFileReader reader = new SAMFileReader(new File(arguments.get(0)));
        reader.setValidationStringency(ValidationStringency.LENIENT);
        SAMRecordIterator it = reader.iterator();

        Bin samBin = new Bin();
        long total = 0;
        try {
          
          while (it.hasNext()) {
            SAMRecord rec = it.next();
            samBin.insert(rec);
            total++;
            if (!rec.getReadUnmappedFlag()) {
              mappedLeft++;
              if (rec.getMappingQuality() >= 1) {
                uniqueLeft++;
                if (rec.getMappingQuality() >= 10) {
                  q10Left++;
                  if (rec.getMappingQuality() >= 20) {
                    q20Left++;
                  }
                }
              }
            }

          }
          
          out("General Resutls for file "  + arguments.get(0) + " :");
          out("Aligner\t%Aligned\t%Unique\t%Q10\t%Q20");
          out(arguments.get(0) + "\t" + mappedLeft / (total) * 100 + "%\t"
              + uniqueLeft / (total) * 100 + "%\t" + q10Left / (total) * 100 + "%\t" + q20Left
              / (total) * 100 + "%");

        } catch (RuntimeEOFException runtimeEOF) {
          err("EOF error!");
          runtimeEOF.printStackTrace();
        }
        out(samBin.toString(total));
      }

      if (findDiscordant) {
        // find discordant reads in each bam file and also show how those discordant reads have been mapped in the other file (alternative mapping)
        if (arguments.size() < 2) {
          out("Please specify at least two files!");
          return;
        }

        File fileRight = new File(arguments.get(1));
        checkFileExists(fileRight);
        
        if (!assumeSorted) {
          SAMUtil.sortByQueryName(fileLeft, TMP_DIR, maxRecordsInRam);
          SAMUtil.sortByQueryName(fileRight, TMP_DIR, maxRecordsInRam);
          fileLeft = new File(arguments.get(0) + ".sorted");
          fileRight = new File(arguments.get(1) + ".sorted");
        }

        SAMFileReader readerLeft = new SAMFileReader(fileLeft);
        SAMFileReader readerLeftCoordinateSorted = null;
        SAMFileReader readerRightCoordinateSorted = null;
        SAMFileReader readerRight = new SAMFileReader(fileRight);

        readerLeft.setValidationStringency(ValidationStringency.LENIENT);
        readerRight.setValidationStringency(ValidationStringency.LENIENT);
        SAMRecordIterator itLeft = readerLeft.iterator();
        SAMRecordIterator itRight = readerRight.iterator();

        SAMRecord rLeft = null, rRight = null;

        String prefix = "";
        
        if (arguments.size() > 2) {
          prefix = arguments.get(2);
        }

        Bin discordantOnlyInLeft =
            new Bin(prefix + fileLeft.getName() + "_DiscordantOnly", true,
                readerLeftCoordinateSorted, outputPair);
        Bin discordantOnlyInRight =
            new Bin(prefix + fileRight.getName() + "_DiscordantOnly", true,
                readerRightCoordinateSorted, outputPair);

        Bin discordantOnlyInLeftOtherPair =
            new Bin(prefix + fileLeft.getName() + "_DiscordantOnly", true,
                readerLeftCoordinateSorted, outputPair);
        Bin discordantOnlyInRightOtherPair =
            new Bin(prefix + fileRight.getName() + "_DiscordantOnly", true,
                readerRightCoordinateSorted, outputPair);

        double logCounterIncrement = 10000;
        double counter = 0;
        double shortReads = 0;
        double longReads = 0;
        double badreads = 0;
        double rightPaired = 0;
        double leftPaired = 0;
        SAMRecord rLeftMate = null, rRightMate = null;
        double leftCorrectlyMapped = 0;
        double rightCorrectlyMapped = 0;

        try {
          while (itLeft.hasNext() && itRight.hasNext()) {
            counter++;
            if (itLeft.hasNext()) {
              rLeft = itLeft.next();
            }
            if (itRight.hasNext()) {
              rRight = itRight.next();
            }
            if (itLeft.hasNext()) {
              rLeftMate = itLeft.next();
            }
            if (itRight.hasNext()) {
              rRightMate = itRight.next();
            }

            if (rLeft.getAlignmentStart() != 0
                && rLeft.getReadName().contains(String.valueOf(rLeft.getAlignmentStart()))) {
              leftCorrectlyMapped++;
            }
            if (rRight.getAlignmentStart() != 0
                && rRight.getReadName().contains(String.valueOf(rRight.getAlignmentStart()))) {
              rightCorrectlyMapped++;
            }

            if ((rLeft.getFirstOfPairFlag() && rRight.getFirstOfPairFlag())
                || (!rLeft.getFirstOfPairFlag() && !rRight.getFirstOfPairFlag())) {

              if (rLeft.getReadLength() < minReadLength) {
                shortReads++;
              } else if (rLeft.getReadLength() > maxReadLength) {
                longReads++;
              } else {
                // discordant bin
                if (Math.abs(rLeft.getAlignmentStart() - rLeft.getMateAlignmentStart()) > maxConcordance
                    && Math.abs(rRight.getAlignmentStart() - rRight.getMateAlignmentStart()) < maxConcordance
                    && !rRight.getReadUnmappedFlag() && !rRight.getMateUnmappedFlag()) {
                  discordantOnlyInLeft.insert(rRightMate);
                  discordantOnlyInLeftOtherPair.insert(rLeftMate);
                  if (rRightMate.getReadName().contains(
                      String.valueOf(rLeftMate.getAlignmentStart()))
                      && !rLeft.getReadUnmappedFlag() && !rLeft.getMateUnmappedFlag()) {
                    rightPaired++;
                    err(rLeft);
                  }
                }

                if (Math.abs(rRight.getAlignmentStart() - rRight.getMateAlignmentStart()) > maxConcordance
                    && Math.abs(rLeft.getAlignmentStart() - rLeft.getMateAlignmentStart()) < maxConcordance
                    && !rLeft.getReadUnmappedFlag() && !rLeft.getMateUnmappedFlag()) {
                  discordantOnlyInRight.insert(rLeftMate);
                  discordantOnlyInRightOtherPair.insert(rRightMate);
                  if (rLeftMate.getReadName().contains(
                      String.valueOf(rRightMate.getAlignmentStart()))) {
                    leftPaired++;
                  }
                }

                if (!rLeft.getReadUnmappedFlag()) {
                  mappedLeft++;
                  if (rLeft.getMappingQuality() >= 1) {
                    uniqueLeft++;
                    if (rLeft.getMappingQuality() >= 10) {
                      q10Left++;
                      if (rLeft.getMappingQuality() >= 20) {
                        q20Left++;
                      }
                    }
                  }
                }

                if (!rRight.getReadUnmappedFlag()) {
                  mappedRight++;
                  if (rRight.getMappingQuality() >= 1) {
                    uniqueRight++;
                    if (rRight.getMappingQuality() >= 10) {
                      q10Right++;
                      if (rRight.getMappingQuality() >= 20) {
                        q20Right++;
                      }
                    }
                  }
                }

              }

              if (counter % logCounterIncrement == 0) {
                out("Processed " + counter + " records so far.");
              }

            } else {
              badreads++;
            }
          }

        } catch (RuntimeEOFException runtimeEOF) {
          err("EOF error!");
          runtimeEOF.printStackTrace();
        }
        discordantOnlyInLeft.close();
        discordantOnlyInRight.close();

        if (outputHTML) {
          DecimalFormat formatter = new DecimalFormat("###.###");

          out(HtmlUtils.htmlHeader());
          out("<table class='percentTable'>");
          out(HtmlUtils.getFiveCellRow("Aligner", "%Aligned", "%Unique", "%Q10", "%Q20"));
          out(HtmlUtils.getFiveCellRow(fileLeft.getName(),
              formatter.format(mappedLeft / (counter - shortReads - longReads) * 100) + "%",
              formatter.format(uniqueLeft / (counter - shortReads - longReads) * 100) + "%",
              formatter.format(q10Left / (counter - shortReads - longReads) * 100) + "%",
              formatter.format(q20Left / (counter - shortReads - longReads) * 100) + "%"));

          out(HtmlUtils.getFiveCellRow(fileRight.getName(),
              formatter.format(mappedRight / (counter - shortReads - longReads) * 100) + "%",
              formatter.format(uniqueRight / (counter - shortReads - longReads) * 100) + "%",
              formatter.format(q10Right / (counter - shortReads - longReads) * 100) + "%",
              formatter.format(q20Right / (counter - shortReads - longReads) * 100) + "%"));

          out("</table>");

          // print output
          String titleLeft = fileLeft.getName();
          String titleRight = fileRight.getName();

          out(HtmlUtils.getBinDiffHtml(discordantOnlyInLeft, discordantOnlyInRight, titleLeft,
              titleRight, true, counter - shortReads - longReads,
              "Discordant Only in One Inverse Stats", false));

          out(HtmlUtils.getBinDiffHtml(discordantOnlyInLeftOtherPair,
              discordantOnlyInRightOtherPair, titleLeft, titleRight, true, counter - shortReads
                  - longReads, "Discordant Only in One Inverse Stats", false));

          out(HtmlUtils.htmlFooter());
          out(fileLeft.getName() + " Paired Wrong: " + rightPaired + " , "
              + fileRight.getName() + " Paired Wrong: " + leftPaired);

          out(fileLeft.getName() + " Correctly Mapped: " + leftCorrectlyMapped
              + " , " + fileRight.getName() + " Correctly Mapped: " + rightCorrectlyMapped);
        }
      }

      if (diff) {
        // compare two BAM files and generate statistics by grouping similar reads into @Bin's
        if (arguments.size() < 2) {
          out("Please specify at least too files!");
          return;
        }

        File fileRight = new File(arguments.get(1));

        if (!assumeSorted) {
          SAMUtil.sortByQueryName(fileLeft, TMP_DIR, maxRecordsInRam);
          SAMUtil.sortByQueryName(fileRight, TMP_DIR, maxRecordsInRam);
          fileLeft = new File(arguments.get(0) + ".sorted");
          fileRight = new File(arguments.get(1) + ".sorted");
        }

        SAMFileReader readerLeft = new SAMFileReader(fileLeft);
        SAMFileReader readerLeftCoordinateSorted = null;
        SAMFileReader readerRightCoordinateSorted = null;

        if (outputPair) {
          readerLeftCoordinateSorted = new SAMFileReader(new File(arguments.get(0) + ".co"));

          readerRightCoordinateSorted = new SAMFileReader(new File(arguments.get(1) + ".co"));
          readerRightCoordinateSorted.setValidationStringency(ValidationStringency.LENIENT);
          readerLeftCoordinateSorted.setValidationStringency(ValidationStringency.LENIENT);

        }
        
        SAMFileReader readerRight = new SAMFileReader(fileRight);

        readerLeft.setValidationStringency(ValidationStringency.LENIENT);
        readerRight.setValidationStringency(ValidationStringency.LENIENT);
        SAMRecordIterator itLeft = readerLeft.iterator();
        SAMRecordIterator itRight = readerRight.iterator();

        SAMRecord rLeft = null, rRight = null;
        Bucketiser boxes = new Bucketiser(maxReadLength, maxReadLength);

        String prefix = "";
        if (arguments.size() > 2) {
          prefix = arguments.get(2);
        }

        Bin sameCigarAndLocationLeft =
            new Bin(prefix + fileLeft.getName() + "_SameCigarAndLocation", true,
                readerLeftCoordinateSorted, false);
        Bin sameCigarAndLocationLeftNoMD =
            new Bin(prefix + fileLeft.getName() + "_SameCigarAndLocationNoMD", true,
                readerLeftCoordinateSorted, false);
        Bin sameCigarAndLocationRightNoMD =
            new Bin(prefix + fileRight.getName() + "_SameCigarAndLocationNoMD", true,
                readerRightCoordinateSorted, outputPair);
        Bin sameCigarAndLocationRight =
            new Bin(prefix + fileRight.getName() + "_SameCigarAndLocation", true,
                readerRightCoordinateSorted, outputPair);
        Bin sameCigarAndFarLocationLeft =
            new Bin(prefix + fileLeft.getName() + "_SameCigarAndFarLocation", true,
                readerLeftCoordinateSorted, false);
        Bin sameCigarAndCloseLocationLeft =
            new Bin(prefix + fileLeft.getName() + "_SameCigarAndCloseLocation", true,
                readerLeftCoordinateSorted, outputPair);
        Bin onlyMappedByRight =
            new Bin(prefix + fileRight.getName() + "_OnlyMapped", true,
                readerRightCoordinateSorted, outputPair);
        Bin differentCigarByLeft =
            new Bin(prefix + fileLeft.getName() + "_DifferentCigar", true,
                readerLeftCoordinateSorted, outputPair);
        Bin sameCigarAndFarLocationRight =
            new Bin(prefix + fileRight.getName() + "_SameCigarAndFarLocation", true,
                readerRightCoordinateSorted, outputPair);
        Bin sameCigarAndCloseLocationRight =
            new Bin(prefix + fileRight.getName() + "_SameCigarAndCloseLocation", true,
                readerRightCoordinateSorted, outputPair);
        Bin onlyMappedByLeft =
            new Bin(prefix + fileLeft.getName() + "_OnlyMapped", true, readerLeftCoordinateSorted,
                outputPair);
        Bin differentCigarByRight =
            new Bin(prefix + fileRight.getName() + "_DifferentCigar", true,
                readerRightCoordinateSorted, outputPair);
        Bin closeLocationDifferentCigarByRight =
            new Bin(prefix + fileRight.getName() + "_CloseLocationDifferentCigar", true,
                readerRightCoordinateSorted, outputPair);
        Bin closeLocationDifferentCigarByLeft =
            new Bin(prefix + fileLeft.getName() + "_CloseLocationDifferentCigar", true,
                readerLeftCoordinateSorted, outputPair);
        Bin sameLocationDifferentCigarByRight =
            new Bin(prefix + fileRight.getName() + "_SameLocationDifferentCigar", true,
                readerRightCoordinateSorted, outputPair);
        Bin sameLocationDifferentCigarByLeft =
            new Bin(prefix + fileLeft.getName() + "_SameLocationDifferentCigar", true,
                readerLeftCoordinateSorted, outputPair);
        Bin unmappedByBoth =
            new Bin(prefix + "UnmappedByBoth", true, readerLeftCoordinateSorted, outputPair);
        Bin discordantOnlyInLeft =
            new Bin(prefix + fileLeft.getName() + "_DiscordantOnly", true,
                readerLeftCoordinateSorted, outputPair);
        Bin discordantOnlyInRight =
            new Bin(prefix + fileRight.getName() + "_DiscordantOnly", true,
                readerRightCoordinateSorted, outputPair);

        double logCounterIncrement = 10000;
        double counter = 0;
        double shortReads = 0;
        double longReads = 0;
        double badreads = 0;
        double sameLocationDifferentCigar = 0;
        double leftDiscordant = 0;
        double rightDiscordant = 0;
        double rightClippings = 0;
        double leftClippings = 0;
        double clippedAndDiscordant = 0;
        double onlyMappedByLeftAndDiscordantInRight = 0;
        double onlyMappedByRightAndDiscordantInLeft = 0;
        double sameCigarAndLocationDifferentMD = 0;
        double leftZeroInMD = 0;
        double rightZeroInMD = 0;


        try {
          while (itLeft.hasNext() && itRight.hasNext()) {
            counter++;
            if (itLeft.hasNext()) {
              rLeft = itLeft.next();
            }
            if (itRight.hasNext()) {
              rRight = itRight.next();
            }

            boxes.insert(rLeft.getReadLength());
            if ((rLeft.getReadLength() != rRight.getReadLength() && (rLeft.getFirstOfPairFlag() == rRight
                .getFirstOfPairFlag()))) {
              err("Same reads in both files do not have the same size!... This might because of sorting issues." + rLeft);
              badreads++;
            }

            if ((rLeft.getFirstOfPairFlag() && rRight.getFirstOfPairFlag())
                || (!rLeft.getFirstOfPairFlag() && !rRight.getFirstOfPairFlag())) {

              if (rLeft.getReadLength() < minReadLength) {
                shortReads++;
              } else if (rLeft.getReadLength() > maxReadLength) {
                longReads++;
              } else {
                // discordant bin
                if (Math.abs(rLeft.getAlignmentStart() - rLeft.getMateAlignmentStart()) > maxConcordance
                    && Math.abs(rRight.getAlignmentStart() - rRight.getMateAlignmentStart()) < maxConcordance
                    && !rRight.getReadUnmappedFlag() && !rRight.getMateUnmappedFlag()) {

                  discordantOnlyInLeft.insert(rLeft);
                }

                if (Math.abs(rRight.getAlignmentStart() - rRight.getMateAlignmentStart()) > maxConcordance
                    && Math.abs(rLeft.getAlignmentStart() - rLeft.getMateAlignmentStart()) < maxConcordance
                    && !rLeft.getReadUnmappedFlag() && !rLeft.getMateUnmappedFlag()) {

                  discordantOnlyInRight.insert(rRight);
                }

                if (!rLeft.getReadUnmappedFlag()) {
                  mappedLeft++;
                  if (rLeft.getMappingQuality() >= 1) {
                    uniqueLeft++;
                    if (rLeft.getMappingQuality() >= 10) {
                      q10Left++;
                      if (rLeft.getMappingQuality() >= 20) {
                        q20Left++;
                      }
                    }
                  }
                }

                if (!rRight.getReadUnmappedFlag()) {
                  mappedRight++;
                  if (rRight.getMappingQuality() >= 1) {
                    uniqueRight++;
                    if (rRight.getMappingQuality() >= 10) {
                      q10Right++;
                      if (rRight.getMappingQuality() >= 20) {
                        q20Right++;
                      }
                    }
                  }
                }

                if (rLeft.getReadName().equals(rRight.getReadName())) {
                  // read exists in both inputs
                  if (rLeft.getReadUnmappedFlag() && rRight.getReadUnmappedFlag()) {
                    unmappedByBoth.insert(rRight);
                  } else if (rLeft.getAlignmentStart() == rRight.getAlignmentStart()
                      && rLeft.getCigar().equals(rRight.getCigar())) {
                    // same cigar AND location
                    sameCigarAndLocationRight.insert(rRight);
                    sameCigarAndLocationLeft.insert(rLeft);
                  } else if (Math.abs(rLeft.getAlignmentStart() - rRight.getAlignmentStart()) < wiggle
                      && rLeft.getCigar().equals(rRight.getCigar())
                      && rLeft.getIntegerAttribute("NM") == rRight.getIntegerAttribute("NM")) {
                    // same cigar and edit distance and CLOSE location
                    sameCigarAndCloseLocationRight.insert(rRight);
                    sameCigarAndCloseLocationLeft.insert(rLeft);
                  } else if (rLeft.getCigar().equals(rRight.getCigar())
                      && rLeft.getIntegerAttribute("NM") == rRight.getIntegerAttribute("NM")) {
                    // same cigar and edit distance FAR location
                    sameCigarAndFarLocationLeft.insert(rLeft);
                    sameCigarAndFarLocationRight.insert(rRight);
                  } else {
                    if (rLeft.getReadUnmappedFlag()) {
                      // not mapped in Left
                      if (rRight.getReadUnmappedFlag()) {
                        // not mapped in either but doen't happen because if both unmapped then
                        // their cigar is the same and this is caught on the above clause but just
                        // in case
                        unmappedByBoth.insert(rRight);
                      } else {
                        // only unmapped by Left
                        onlyMappedByRight.insert(rRight);
                        if (!rLeft.getMateUnmappedFlag()) {
                          onlyMappedByRightAndDiscordantInLeft++;
                        }
                      }
                    } else {
                      if (rRight.getReadUnmappedFlag()) {
                        // unmapped by Right only
                        onlyMappedByLeft.insert(rLeft);
                        if (!rRight.getMateUnmappedFlag()) {
                          onlyMappedByLeftAndDiscordantInRight++;
                        }
                      } else {
                        //same location but different cigar
                        if (rLeft.getAlignmentStart() == rRight.getAlignmentStart()) {
                          sameLocationDifferentCigarByLeft.insert(rLeft);
                          sameLocationDifferentCigarByRight.insert(rRight);
                        } else if (Math.abs(rLeft.getAlignmentStart() - rRight.getAlignmentStart()) < 30) {
                          //close location different cigar
                          closeLocationDifferentCigarByLeft.insert(rLeft);
                          closeLocationDifferentCigarByRight.insert(rRight);
                        } else {
                          // different cigar and mapped by both
                          differentCigarByLeft.insert(rLeft);
                          differentCigarByRight.insert(rRight);
                        }
                      }
                    }
                  }
                } else {
                  // read does not exist in both samples!!!
                  err(rLeft.getReferenceName() + "," + rLeft.getReadName()
                      + " Files do not match.");
                }
              }

              if (counter % logCounterIncrement == 0) {
                System.out.println("Processed " + counter + " records so far.");
              }

            } else {
              badreads++;
            }
          }

        } catch (RuntimeEOFException runtimeEOF) {
          err("EOF error!");
          runtimeEOF.printStackTrace();
        }

        // flush bins
        sameCigarAndCloseLocationLeft.close();
        sameCigarAndCloseLocationRight.close();
        sameCigarAndFarLocationLeft.close();
        sameCigarAndFarLocationRight.close();
        sameCigarAndLocationLeft.close();
        sameCigarAndLocationRight.close();
        differentCigarByLeft.close();
        differentCigarByRight.close();
        unmappedByBoth.close();
        discordantOnlyInLeft.close();
        discordantOnlyInRight.close();
        onlyMappedByRight.close();
        onlyMappedByLeft.close();

        if (outputHTML) {
          //generate HTML report
          DecimalFormat formatter = new DecimalFormat("###.###");

          out(HtmlUtils.htmlHeader());
          out("Same Cigar and Location Different MD that have N's:"
              + sameCigarAndLocationDifferentMD);

          out("Only mapped in " + fileRight.getName() + " and discordant in "
              + fileLeft.getName() + "  : " + onlyMappedByRightAndDiscordantInLeft);
          out("Only mapped in " + fileLeft.getName() + " and discordant in "
              + fileRight.getName() + "  : " + onlyMappedByLeftAndDiscordantInRight);

          out("<table class='percentTable'>");
          out(HtmlUtils.getFiveCellRow("Aligner", "%Aligned", "%Unique", "%Q10", "%Q20"));
          out(HtmlUtils.getFiveCellRow(fileLeft.getName(),
              formatter.format(mappedLeft / (counter - shortReads - longReads) * 100) + "%",
              formatter.format(uniqueLeft / (counter - shortReads - longReads) * 100) + "%",
              formatter.format(q10Left / (counter - shortReads - longReads) * 100) + "%",
              formatter.format(q20Left / (counter - shortReads - longReads) * 100) + "%"));

          out(HtmlUtils.getFiveCellRow(fileRight.getName(),
              formatter.format(mappedRight / (counter - shortReads - longReads) * 100) + "%",
              formatter.format(uniqueRight / (counter - shortReads - longReads) * 100) + "%",
              formatter.format(q10Right / (counter - shortReads - longReads) * 100) + "%",
              formatter.format(q20Right / (counter - shortReads - longReads) * 100) + "%"));

          out("</table>");

          // print output
          String titleLeft = fileLeft.getName();
          String titleRight = fileRight.getName();

          out(HtmlUtils.getBinDiffHtml(sameCigarAndLocationLeft, sameCigarAndLocationRight,
              titleLeft, titleRight, true, counter - shortReads - longReads,
              "Same Cigar And Location", false));
          out(HtmlUtils.getBinDiffHtml(sameCigarAndFarLocationLeft,
              sameCigarAndFarLocationRight, titleLeft, titleRight, true, counter - shortReads
                  - longReads, "Same Cigar And Far Location", false));
          out(HtmlUtils.getBinDiffHtml(sameCigarAndCloseLocationLeft,
              sameCigarAndCloseLocationRight, titleLeft, titleRight, true, counter - shortReads
                  - longReads, "Same Cigar And Close Location", false));

          out(HtmlUtils.getBinDiffHtml(onlyMappedByLeft, onlyMappedByRight, titleLeft,
              titleRight, true, counter - shortReads - longReads, "Only Mapped by One", false));

          out(HtmlUtils.getBinDiffHtml(differentCigarByLeft, differentCigarByRight, titleLeft,
              titleRight, true, counter - shortReads - longReads,
              "Different Cigar Different Location", false));
          out(HtmlUtils.getBinDiffHtml(closeLocationDifferentCigarByLeft,
              closeLocationDifferentCigarByRight, titleLeft, titleRight, true, counter - shortReads
                  - longReads, "Different Cigar Close Location", false));

          out(HtmlUtils.getBinDiffHtml(sameLocationDifferentCigarByLeft,
              sameLocationDifferentCigarByRight, titleLeft, titleRight, true, counter - shortReads
                  - longReads, "Same Location Different Cigar", false));

          out(HtmlUtils.getBinDiffHtml(unmappedByBoth, unmappedByBoth, titleLeft, titleRight,
              true, counter - shortReads - longReads, "Unmapped By Both", false));

          out(HtmlUtils.getBinDiffHtml(discordantOnlyInLeft, discordantOnlyInRight, titleLeft,
              titleRight, true, counter - shortReads - longReads, "Discordant Only in One", false));

          out(HtmlUtils.getBinDiffHtml(sameCigarAndLocationLeftNoMD,
              sameCigarAndLocationRightNoMD, titleLeft, titleRight, true, counter - shortReads
                  - longReads, "Same CIGAR and Location different MD", false));

          out("For Different CIGAR and Far Location, discordant reads in "
              + fileLeft + ": " + leftDiscordant + "<br>");
          out("Discordant reads in " + fileRight + ": " + rightDiscordant + "<br>");
          out("Clipped reads in " + fileLeft + ": " + leftClippings + "<br>");
          out("Clipped reads in " + fileRight + ": " + rightClippings + "<br>");
          out("Clipped and discordant in " + fileRight + " : "
              + clippedAndDiscordant + "<br>");
          out("Bad reads:" + badreads + "(" + badreads / counter * 100 + "%)"
              + "<br>");

          out("Short reads:" + shortReads + "(" + shortReads / counter * 100 + "%)"
              + "<br>");

          out("Long reads:" + longReads + "(" + longReads / counter * 100 + "%)"
              + "<br>");

          out("Same location different CIGAR" + sameLocationDifferentCigar + "<br>");

          out("Only considering reads longer than " + minReadLength
              + " and shorted than " + maxReadLength + " inclusive." + "<br>");
          out(fileLeft.getName() + " Zero in MD" + leftZeroInMD);
          out(fileRight.getName() + " Zero in MD" + rightZeroInMD);
          out(HtmlUtils.htmlFooter());

        } else {
          out("Read distribution: " + boxes);
          out("Bad reads:" + badreads + "(" + badreads / counter * 100 + "%)");

          out("Short reads:" + shortReads + "(" + shortReads / counter * 100 + "%)");

          out("Long reads:" + longReads + "(" + longReads / counter * 100 + "%)");

          out("Same location different CIGAR" + sameLocationDifferentCigar);

          out("Only considering reads longer than " + minReadLength
              + " and shorted than " + maxReadLength + " inclusive.");

          DecimalFormat formatter = new DecimalFormat("###.###");

          out("General Results:");
          out("Aligner\t%Aligned\t%Unique\t%Q10\t%Q20" );
          out(fileLeft.getName() + "\t"
                  + formatter.format(mappedLeft / (counter - shortReads - longReads) * 100) + "%\t"
                  + formatter.format(uniqueLeft / (counter - shortReads - longReads) * 100) + "%\t"
                  + formatter.format(q10Left / (counter - shortReads - longReads) * 100) + "%\t"
                  + formatter.format(q20Left / (counter - shortReads - longReads) * 100) + "%");
          out(fileRight.getName() + "\t"
              + formatter.format(mappedRight / (counter - shortReads - longReads) * 100) + "%\t"
              + formatter.format(uniqueRight / (counter - shortReads - longReads) * 100) + "%\t"
              + formatter.format(q10Right / (counter - shortReads - longReads) * 100) + "%\t"
              + formatter.format(q20Right / (counter - shortReads - longReads) * 100) + "%");
          
          out(sameCigarAndLocationLeft.toString(counter - shortReads - longReads));
          out(sameCigarAndLocationRight.toString(counter - shortReads - longReads));

          out(sameCigarAndFarLocationLeft.toString(counter - shortReads - longReads));
          out(sameCigarAndFarLocationRight
              .toString(counter - shortReads - longReads));

          out(sameCigarAndCloseLocationLeft.toString(counter - shortReads
              - longReads));
          out(sameCigarAndCloseLocationRight.toString(counter - shortReads
              - longReads));

          out(onlyMappedByLeft.toString(counter - shortReads - longReads));
          out(onlyMappedByRight.toString(counter - shortReads - longReads));

          out(differentCigarByLeft.toString(counter - shortReads - longReads));
          out(differentCigarByRight.toString(counter - shortReads - longReads));

          out(closeLocationDifferentCigarByLeft.toString(counter - shortReads
              - longReads));
          out(closeLocationDifferentCigarByRight.toString(counter - shortReads
              - longReads));

          out(sameLocationDifferentCigarByLeft.toString(counter - shortReads
              - longReads));
          out(sameLocationDifferentCigarByRight.toString(counter - shortReads
              - longReads));

          out(unmappedByBoth.toString(counter - shortReads - longReads));

          out("For Different CIGAR and Far Location, discordant reads in "
              + fileLeft + ": " + leftDiscordant);
          out("Discordant reads in " + fileRight + ": " + rightDiscordant);
          out("Clipped reads in " + fileLeft + ": " + leftClippings);
          out("Clipped reads in " + fileRight + ": " + rightClippings);
          out("Clipped and discordant in " + fileRight + " : "
              + clippedAndDiscordant);

        }
      }

    } catch (CmdLineException e) {
      parser.printUsage(System.err);
      return;
    }
  }
  
  /** Wrapper to save typing :) */
  private void out(Object obj){
    System.out.println(obj);
  }
  
  /** Wrapper to save typing :) */
  private void err(Object obj){
    throw new RuntimeException(obj.toString());
  }

  private void checkFileExists(File file){
    if (!file.exists()){
      err("File " + file + " does not exist!");
    }
  }
}
