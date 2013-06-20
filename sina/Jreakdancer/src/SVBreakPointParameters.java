import java.util.HashMap;


public class SVBreakPointParameters {

  public String sv_chr1 = "", sv_chr2 = "";
  public Integer    sv_pos1 = 0, sv_pos2 = 0;
  public String sv_ori1, sv_ori2;
  
  public int normal_rp;

   public HashMap<String, Integer> read_count=new HashMap<String, Integer>();
  
  public Integer firstNode=0;
  
  public double copy_number;
  public String SVType="UN";
  public Integer PhredQ=0;
  
}
