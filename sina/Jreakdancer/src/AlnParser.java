import java.awt.image.SampleModel;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import net.sf.picard.sam.SamPairUtil;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMUtils;

public class AlnParser {

  public static void main(String[] args) {
    
    String inputFile = "/Users/sina/tmp/chr1.bam";
    SAMFileReader r = new SAMFileReader(new File(inputFile));
    int count = 0;
    SAMRecordIterator iter = r.iterator();
    while (count++ < 10000000) {
      SAMRecord record = iter.next();
      System.out.println(record.getReadName());
      // alnParse(record);
    }
    iter.close();
    r.close();
  }

  public static int alnParse(SAMRecord record) {
    if (record.getDuplicateReadFlag()) return 0;
    String ori = "+";
    if (record.getReadNegativeStrandFlag()) ori = "-";
    if (record.getReadPairedFlag()) {
      if (record.getReadUnmappedFlag()) return 192;
      if (record.getMateUnmappedFlag()) return 64;
      if (!record.getReferenceName().equals(record.getMateReferenceName())) return 32;
      String ori2 = "+";
      if (record.getMateNegativeStrandFlag()) ori2 = "-";
      if (record.getProperPairFlag()) {
        if (record.getAlignmentStart() < record.getMateAlignmentStart()) {
          if (ori.equals("+")) return 18;
        } else {
          if (ori.equals("-")) return 18;
        }
        return 20;
      } else {
        if (ori.equals(ori2)) {
          if (ori.equals("+")) return 1;
          return 8;
        }
        if (record.getMateAlignmentStart() > record.getAlignmentStart()) {
          if (ori.equals("-")) return 4;
          return 2;
        } else {
          if (ori.equals("+")) return 4;
          return 2;
        }

      }
    }
    return 0;
  }
}
