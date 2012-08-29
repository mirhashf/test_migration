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

import java.text.DecimalFormat;

/**
 * Very simple HTML methods to generate the report
 * 
 * @author amir@binatechnologies.com
 *
 */
public class HtmlUtils {
  public static String getFiveCellRow(String c1, String c2, String c3, String c4, String c5) {
    String out = "";
    out +=
        "<tr><td>" + c1 + "</td><td>" + c2 + "</td><td>" + c3 + "</td><td>" + c4 + "</td><td>" + c5
            + "</td></tr>";
    return out;
  }

  public static String getBinDiffHtml(Bin leftBin, Bin rightBin, String titleLeft, String titleRight,
      boolean printSubBin, double total, String binTitle, boolean isUnmappedBin) {
    String out = "";
    if (printSubBin) {
      out += "<table class='diffTable'>";
      out += getSuperBinRow(binTitle);
      out += getTitleRowHtml(titleLeft, titleRight);
      out += getRowHtmlWithPercent("Total Bin Size", leftBin.total, rightBin.total, total, total);
      out += getBinRow("MAPQ>5 Sub-Bin");

    }

    out += getRowHtmlWithPercent("Bin size", leftBin.CountReads, rightBin.CountReads, total, total);

    if (!isUnmappedBin) {
      out +=
          getRowHtmlWithAverage("Sum Gap Open", leftBin.SumGapOpen, rightBin.SumGapOpen,
              leftBin.CountReads - leftBin.CountUnmapped, rightBin.CountReads
                  - rightBin.CountUnmapped);
      out +=
          getRowHtmlWithAverage("Sum Edit Distance", leftBin.SumEditDistance,
              rightBin.SumEditDistance, leftBin.CountReads - leftBin.CountUnmapped,
              rightBin.CountReads - rightBin.CountUnmapped);
      out +=
          getRowHtmlWithAverage("Sum Mismatches", leftBin.SumMismatches, rightBin.SumMismatches,
              leftBin.CountReads - leftBin.CountUnmapped, rightBin.CountReads
                  - rightBin.CountUnmapped);
      out +=
          getRowHtmlWithAverage("Sum Gap Extensions", leftBin.SumGapExtensions,
              rightBin.SumGapExtensions, leftBin.CountReads - leftBin.CountUnmapped,
              rightBin.CountReads - rightBin.CountUnmapped);
      out +=
          getRowHtmlWithAverage("Sum Mapping Quality", leftBin.SumMappingQuality,
              rightBin.SumMappingQuality, leftBin.CountReads - leftBin.CountUnmapped,
              rightBin.CountReads - rightBin.CountUnmapped);
      out +=
          getRowHtmlWithAverage("Sum Soft Clipping", leftBin.SumSoftClipping,
              rightBin.SumSoftClipping, leftBin.CountReads - leftBin.CountUnmapped,
              rightBin.CountReads - rightBin.CountUnmapped);
      out +=
          getRowHtmlWithAverage("Sum Normalized Edit Distance", leftBin.SumNormalizedEditDistance,
              rightBin.SumNormalizedEditDistance, leftBin.CountReads - leftBin.CountUnmapped,
              rightBin.CountReads - rightBin.CountUnmapped);
      out +=
          getRowHtmlWithAverage("Sum of 4-mers Repetitiveness", leftBin.SumMostRepetetiveKmers,
              rightBin.SumMostRepetetiveKmers, leftBin.CountReads, rightBin.CountReads);
    }

    out +=
        getRowHtmlWithAverage("Sum Base Quality", leftBin.SumBaseQuality, rightBin.SumBaseQuality,
            leftBin.CountBases, rightBin.CountBases);
    out +=
        getRowHtmlWithAverage("Sum Read Length", leftBin.CountBases, rightBin.CountBases,
            leftBin.CountReads - leftBin.CountUnmapped, rightBin.CountReads
                - rightBin.CountUnmapped);
    if (!isUnmappedBin) {
      out += getRowHtml("Number of Discordant", leftBin.discordantReads, rightBin.discordantReads);
    }
    if (printSubBin) {
      out += getBinRow("MAPQ<5 Sub-Bin");
      out +=
          getBinDiffHtml(leftBin.lowMapqBin, rightBin.lowMapqBin, "", "", false, total, "", false);

      out += getBinRow("Unmapped Sub-Bin");
      out +=
          getBinDiffHtml(leftBin.unmappedReadsBin, rightBin.unmappedReadsBin, "", "", false, total,
              "", true);
    }

    if (printSubBin) {
      out += "</table>";
    }

    return out;
  }

  public static String getBinRow(String binName) {
    return "<tr><td colspan='3' class='binTitle'>" + binName + "</td></tr>";
  }

  public static String getSuperBinRow(String binName) {
    return "<tr><td colspan='3' class='superBinTitle'>" + binName + "</td></tr>";
  }

  public static String htmlHeader() {
    String out = "<html><header><title>SeqAlto Alignment Comparator</title>";
    out +=
        "<style>body,td {font-family: 'Helvetica Neue';font-size: 12px;} "
            + ".leftVal{background:#EEE} "
            + ".superBinTitle{color: #DD4B39;"
            + "font-size: 12px;border-bottom: 1px solid #CCC;"
            + "height: 30px;font-weight: bold;text-align: center;} "
            + ".rowTitle {backgroun:white} "
            + ".diffTable td{ width:200px;padding: 3px;} "
            + ".diffTable{margin-bottom: 10px; border: 1px solid #AAA} "
            + ".binTitle{background: #999;color: white;text-align: center;} "
            + ".chart div {   font: 10px sans-serif;   "
            + "background-color: steelblue;   "
            + "text-align: right;   padding: 3px;   margin: 1px;   "
            + "color: white;   display: inline-block; } .percentTable {width:628px; border: solid #CCC 1px} .percentTable td {border: solid #CCC 1px} "
            + ".bar{   font: 10px sans-serif;    display: block; }" + "</style></header><body>";
    return out;
  }

  public static String htmlFooter() {
    String out = "</body>";
    return out;
  }

  public static String getRowHtmlWithAverage(String title, double left, double right, double leftTotal,
      double rightTotal) {

    DecimalFormat formatter = new DecimalFormat("###.###");
    String out = "";
    out +=
        getRowHtml(title, formatter.format(left) + " (ave: " + formatter.format(left / leftTotal)
            + ")", formatter.format(right) + " (ave: " + formatter.format(right / rightTotal) + ")");
    return out;
  }

  public static String getRowHtmlWithPercent(String title, double left, double right, double leftTotal,
      double rightTotal) {
    DecimalFormat formatter = new DecimalFormat("###.###");
    String out = "";
    out +=
        getRowHtml(title, left + " (" + formatter.format(left / leftTotal * 100) + "%)", right
            + " (" + formatter.format(right / rightTotal * 100) + "%)");
    return out;
  }

  public static String getRowHtml(String title, Object left, Object right) {
    String out = "";
    out +=
        "<tr><td class='rowTitle'>" + title + "</td><td class='leftVal'>" + left
            + "</td><td class='rightVal'>" + right + "</td></tr>";
    return out;
  }

  public static String getTitleRowHtml(Object left, Object right) {
    String out = "";
    out +=
        "<tr><td></td><td class='leftVal'>" + left + "</td><td class='rightVal'>" + right
            + "</td></tr>";
    return out;
  }
}
