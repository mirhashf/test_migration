import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;


public class FRecord {
  
  public FRecord(SAMRecord record){
    this.record=record;
  }
  
  public SAMRecord record;
  public int flagZ=0;
}
