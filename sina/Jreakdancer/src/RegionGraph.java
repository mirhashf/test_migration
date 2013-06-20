import java.util.HashMap;


public class RegionGraph extends HashMap<Integer,HashMap<Integer,Integer>> {

  public void add(Integer region, Integer region2){
    add(region,region2,1);
  }

  public void add(Integer region, Integer region2, int numLinks){
    addAsym(region,region2,numLinks);
    addAsym(region2,region,numLinks);
  }
  private void addAsym(Integer region, Integer region2, int numLinks){
    if(!containsKey(region))
      put(region,new HashMap<Integer,Integer>());
    HashMap<Integer,Integer> secondNode=this.get(region);
    if(!secondNode.containsKey(region2))
      secondNode.put(region2, 0);
    secondNode.put(region2,secondNode.get(region2)+numLinks);
  }
  
  public void remove(Integer region,Integer region2){
    removeAsym(region,region2);
    removeAsym(region2,region);
  }

  private void removeAsym(Integer region, Integer region2) {
    if(!containsKey(region))
      return;
    HashMap<Integer, Integer> secondNode = get(region);
    if(!secondNode.containsKey(region2))
      return;
    secondNode.put(region2,0);
  }
  

}