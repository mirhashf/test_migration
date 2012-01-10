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

import java.text.DecimalFormat;

/**
 * Puts input into intervals (buckets) to be used for statistics. (horrible naming probably :D )
 * 
 * @author amir@binatechnologies.com (Amirhossein Kiani)
 */
public class Bucketiser {
  int numBuckets = 100;
  long[] buckets;
  double range = 0;
  double total = 0;
  // to format numbers
  DecimalFormat formatter = new DecimalFormat("###.###");

  /* Creates a new @Bucketiser instance */
  public Bucketiser(double range, int numIntervals) {
    buckets = new long[numIntervals + 1];
    numBuckets = numIntervals;
    this.range = range;
  }

  /* Inserts a new number to the @Bucketiser putting it in the right bucket */
  public void insert(double number) {
    buckets[(int) (number / (range / numBuckets))]++;
    total++;
  }

  /* Writes out the contents */
  @Override
  public String toString() {
    String out = "";
    for (int i = 0; i < buckets.length; i++) {

      out += i * (range / numBuckets) + "\t " + buckets[i] + "\n";
    }
    return out.toString();
  }

  /* Convert contents to JSON (used for generating JS plots) */
  public String getJsArray() {
    String out = "[";
    for (int i = 0; i < buckets.length; i++) {
      out += formatter.format(buckets[i] / total * 100) + ",";
    }
    out += "]";
    return out.toString();
  }
}
