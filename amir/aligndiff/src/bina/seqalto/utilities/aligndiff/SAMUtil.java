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

import java.io.File;

import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriterImpl;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
/**
 * General methods to deal with SAM/BAM data
 * 
 * @author amir@binatechnologies.com (Amirhossein Kiani)
 *
 */
public class SAMUtil {
  public static int getIntegerFlagSafe(SAMRecord samRecord, String flag) {

    if (samRecord.getIntegerAttribute(flag) != null)
      return samRecord.getIntegerAttribute(flag);
    else {
      System.err.println("No " + flag + " for read " + samRecord.getReadName());
      return 0;
    }
  }

  /* Sorts the BAM/SAM file by query name(Read ID) */
  public static void sortByQueryName(File file, File TMP_DIR, int maxRecordsInRam) {
    BAMFileWriter writer = new BAMFileWriter(new File(file.getPath() + ".sorted"));
    SAMFileWriterImpl.setDefaultMaxRecordsInRam(maxRecordsInRam);

    if (!TMP_DIR.exists()) {
      // Intentionally not checking the return value, because it may be that
      // the program does not
      // need a tmp_dir. If this fails, the problem will be discovered
      // downstream.
      TMP_DIR.mkdirs();
    }
    System.setProperty("java.io.tmpdir", TMP_DIR.getAbsolutePath());

    writer.setSortOrder(SAMFileHeader.SortOrder.queryname, false);
    SAMFileHeader header = new SAMFileHeader();
    SAMFileReader reader = new SAMFileReader(file);
    reader.setValidationStringency(ValidationStringency.LENIENT);
    header.setSequenceDictionary(reader.getFileHeader().getSequenceDictionary());
    writer.setHeader(header);
    SAMRecordIterator it = reader.iterator();
    while (it.hasNext()) {
      writer.addAlignment(it.next());
    }
    writer.close();
  }

}
