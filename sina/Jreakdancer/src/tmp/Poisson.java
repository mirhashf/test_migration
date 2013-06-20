package tmp;

public class Poisson {

  public static void main(String[] args){
    int n=0;
      System.out.println(n);
  }
  
   
  private double PI = Math.PI;
  private double LZERO = -1e10;
  private double LSMALL = LZERO/2;
  private double SMALL = Math.exp(LSMALL);
  private double minLogExp = -Math.log(-LZERO);

  protected double LogPoissonTailProb(float n, float lambda){
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
