package bina.seqalto.utilities.bamextractor;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

public class ReadBamFile {
  public static void main(String[] args) {
    String bamPath = "/Users/sina/git/seqalto/variant-caller/reduoso/test/resources/bam/ds3.seqalto.short.0.bam";
    
    SAMFileReader testReader = new SAMFileReader(new File(bamPath));
    for (SAMRecord rec : testReader){
      System.out.println(rec.getReadString());
    }
    testReader.close();
  }
}
