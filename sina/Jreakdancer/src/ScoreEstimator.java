import java.util.ArrayList;


public class ScoreEstimator {

  public SVDetector svDetector;

  public ScoreEstimator(SVDetector svDetector) {
    this.svDetector = svDetector;
  }

  

  private int totalBadRegionSize(ArrayList<Integer> regions) {
    int total_region_size = 0;
    for (Integer region : regions) {

      int clust_start = svDetector.an.ap.regionInfo.regs_name.get(region).beginc;
      int clust_end = svDetector.an.ap.regionInfo.regs_name.get(region).lastc;
      total_region_size += clust_end - clust_start + 1;
    }
    return total_region_size;
  }



  protected double computeProbScore(ArrayList<Integer> regions, int flag) {
   try{
    double refLen = svDetector.an.cp.refLen;
    int total_region_size = totalBadRegionSize(regions);

    double lambda;
    double logpvalue = 0;
    for (String rg : svDetector.an.cp.numAllBadReadsPerType.get(flag).keySet()) {
      lambda =
          ((double) total_region_size) * svDetector.an.cp.numAllBadReadsPerType.get(flag).get(rg)
              / refLen;
     logpvalue += LogPoissonTailProb((double)svDetector.SVTypeLibraryReadcount.get(flag).get(rg), lambda);
    }
    return logpvalue;
   }catch(Exception e){ return 0;}
  }

  private double PI = Math.PI;
  private double LZERO = -1e10;
  private double LSMALL = LZERO/2;
  private double SMALL = Math.exp(LSMALL);
  private double minLogExp = -Math.log(-LZERO);

  protected double LogPoissonTailProb(double n, double lambda){
      double logprob = LZERO;
      double plogprob;
      do{
          plogprob = logprob;
          logprob = LAdd(plogprob, LogPoissonPDF(n++,lambda));
      } while(logprob-plogprob < 0.01);
      return logprob;
  }

  protected double LogPoissonPDF(double k, double lambda){
      double logk_factorial = k==0 ? 0 : k * (Math.log (k) ) - k + 0.5 * (Math.log (2*PI*k) ) ;
      double log_lambda = lambda<=0 ? LSMALL:Math.log(lambda);
      double logp = k*log_lambda - lambda - logk_factorial;
      return logp;
  }

  protected double LAdd(double x, double y){
      double temp, diff, z;
      if(x < y){
          temp = x;
          x = y;
          y = temp;
      }
      diff = y - x;
      if(diff < minLogExp)
          return x<LSMALL?LZERO:x;
      else{
          z = Math.exp(diff);
          return x+Math.log(1.0+z);
      }
     
  }

}
