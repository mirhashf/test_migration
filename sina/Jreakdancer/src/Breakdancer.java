import java.io.File;
import java.io.IOException;


public class Breakdancer {

  public static int MINMAPQ = 35;
  public static int maxEventSize= 1000000000;
  
  
  public static void main(String[] args) throws IOException {
    args = new String[] {"/Users/sina/tmp/chr1.bam", "/Users/sina/tmp/chr1.stats","/Users/sina/tmp/scratch", "chr1","0","0","3",MINMAPQ+""};
     String output = args[1];
     if (new File(output).exists()==false) 
      BamStatistics.getBamStatistics(args);
    
     String scratchDir=args[2];
     if(new File(scratchDir).list()==null || new File(scratchDir).list().length==0)
       Candidator.getBamStatistics(args);
     
     Analyzer.getAnalyzer(args);
     
    
     
  }
}
