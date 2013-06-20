package tmp;
public class LibInfo {

  public LibInfo(String map, double mean, double std, double upper, double lower, int readLen) {
    super();
    this.map = map;
    this.mean = mean;
    this.std = std;
    this.upper = upper;
    this.lower = lower;
    this.readLen = readLen;
  }

  @Override
  public String toString() {
    return "LibInfo [map=" + map + ", mean=" + mean + ", std=" + std + ", upper=" + upper
        + ", lower=" + lower + ", readLen=" + readLen + "]";
  }

  String map;
  double mean;
  double std;
  double upper;
  double lower;
  int readLen;
}
