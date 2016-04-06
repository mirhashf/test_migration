/**
 * Roman Number Adder
 * @author Kuang Sheng
 * @summary this program add two Roman numerals (in String) without converting
 * 			them to decimal or any base.
 */

import java.util.*;

public class romanNumberAdder {
	
	private final String order = "IVXLCDM";
	private final String[] subtractive = new String[]{
		"IV"  ,"IX"   ,"XL"  ,"XC"   ,"CD"  ,"CM"};
	private final String[] additive = new String[]{
		"IIII","VIIII","XXXX","LXXXX","CCCC","DCCCC"};
	private final String[] additivePlus = new String[]{
		"IIIII","VV","XXXXX","LL","CCCCC","DD"};
	private final String[] subtractivePlus = new String[]{
		"V"    ,"X" ,"L"    ,"C" ,"D"    ,"M"};
	
	//used for conversion from subtractive to additive form
	private HashMap<String,String> convertMap;
	//used for conversion from additive to combined form
	private HashMap<String,String> convertBackMap;
	
	public romanNumberAdder(){
		convertMap = new HashMap<String,String>();
		for(int i=0; i<subtractive.length; i++){
			convertMap.put(subtractive[i],additive[i]);
		}
		convertBackMap = new HashMap<String,String>();
		for(int i=0; i<additivePlus.length; i++){
			convertBackMap.put(additivePlus[i],subtractivePlus[i]);
		}
	}
	
	public static void main(String[] args) {
		String num1 = "";
		String num2 = "";
		if(args!=null && args.length>=2){
			num1 = args[0];
			num2 = args[1];
		}
		romanNumberAdder rna = new romanNumberAdder();
		//System.out.println(rna.add("MCLXXIV","CXXXIX"));
		System.out.println(rna.add(num1,num2));
	}
	
	/**
	 * add()
	 * @param r1: Roman number string 1
	 * @param r2: Roman number string 2
	 * @return added result Roman number in String
	 */
	public String add(String r1, String r2){
		
		//special input cases: null empty
		if(r1==null && r2==null) return "";
		if(r1==null || r2==null) return r1==null?r2:r1;
		if(r1.length()==0 || r2.length()==0) return r1.length()==0?r2:r1;
		
		r1 = r1.toUpperCase();
		r2 = r2.toUpperCase();

		String r1_additive = convertSubtractiveToAdditive(r1);
		String r2_additive = convertSubtractiveToAdditive(r2);
		String combineR1R2str = r1_additive+r2_additive;
		
		//need Character array to work with myComparator
		Character[] combineR1R2 = new Character[r1_additive.length()+r2_additive.length()];
		for(int i=0; i<combineR1R2str.length(); i++){
			combineR1R2[i]=combineR1R2str.charAt(i);
		}
		
		Arrays.sort(combineR1R2,new myComparator());
		
		String res = combineAdditive(characterArrayToString(combineR1R2));
		
		return res;
	}
	
	/**
	 * myComparator
	 * feed in Arrays.sort, sort the Character in descending order based on 
	 * the order constant string defined above.
	 *
	 */
	class myComparator implements Comparator<Character>{
		@Override
		public int compare(Character o1, Character o2) {
			return order.indexOf(Character.toUpperCase(o2))
				   -order.indexOf(Character.toUpperCase(o1));
		}
	}
	
	/**
	 * convertSubtractiveToAdditive()
	 * @param num: Roman number string
	 * @return modified Roman number string 
	 * There are 6 cases of subtractive forms:
	 * 	"IV"  ,"IX"   ,"XL"  ,"XC"   ,"CD"  ,"CM"
	 * convert to additive forms:
	 * 	"IIII","VIIII","XXXX","LXXXX","CCCC","DCCCC"
	 */
	public String convertSubtractiveToAdditive(String num){
		int i=0;
		while(i<num.length()-1){
			char cur = num.charAt(i);
			char next = num.charAt(i+1);
			if(order.indexOf(cur)<order.indexOf(next)){
				String toAdditive = convertMap.get(""+cur+next);
				num = num.substring(0,i)+toAdditive+num.substring(i+2,num.length());
			}
			i++;
		}
		return num;
	}
	
	/**
	 * combineAdditive()
	 * @param num: Roman number string
	 * @return modified Roman number string 
	 * There are 6 cases of additive forms:
	 * 	"IIIII","VV","XXXXX","LL","CCCCC","DD"
	 * convert to combined forms:
	 * 	"V"    ,"X" ,"L"    ,"C" ,"D"    ,"M"
	 */
	public String combineAdditive(String num){
		for(int i=0; i<additivePlus.length; i++){
			String toChange = additivePlus[i];
			int pos = num.indexOf(toChange);
			while(pos>=0){
				num = num.substring(0,pos)+convertBackMap.get(toChange)
					  +num.substring(pos+toChange.length(),num.length());
				pos = num.indexOf(additivePlus[i]);
			}
		}
		return num;
	}
	
	/**
	 * characterArrayToString()
	 * @param strs: array of Characters
	 * @return string concatenated by elements in array
	 */
	public String characterArrayToString(Character[] strs){
		StringBuilder sb = new StringBuilder();
		for(int i=0; i<strs.length; i++){
			sb.append(strs[i]);
		}
		return sb.toString();
	}

}
