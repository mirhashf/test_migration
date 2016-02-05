import java.util.Arrays;
import java.util.Comparator;
import java.util.Scanner;

public class RomanAddition {
	
	//Substitution before adding
	public static String substitute(String s){
		s = s.replaceAll("IV", "IIII");
		s = s.replaceAll("IX", "VIIII");
	    s = s.replaceAll("XL", "XXXX");
	    s = s.replaceAll("XC", "LXXXX");
	    s = s.replaceAll("CD", "CCCC");
	    s = s.replaceAll("CM", "DCCCC");
	    return s;
	}
	
	//Substitution after adding
	public static String reverse_substitute(String s){
		s = s.replaceAll("IIIII", "V");
		s = s.replaceAll("VV", "X");
		s = s.replaceAll("XXXXX", "L");
		s = s.replaceAll("LL", "C");
		s = s.replaceAll("CCCCC", "D");
		s = s.replaceAll("DD", "M");
		s = s.replaceAll("VIIII", "IX");
		s = s.replaceAll("IIII", "IV");
		s = s.replaceAll("LXXXX", "XC");
		s = s.replaceAll("XXXX", "XL");
		s = s.replaceAll("DCCCC", "CM");
		s = s.replaceAll("CCCC", "CD");
		
		return s;
	}
	
	//Sorting as per roman format
	public static String sort(String sum){
		Character[] sum1 = new Character[sum.length()];
		String custom = "MDCLXVI";
		
		for(int i=0;i<sum1.length;i++){
			sum1[i] = sum.charAt(i);
		}
		
		Arrays.sort(sum1, new Comparator<Character>(){
			@Override
			public int compare(Character s1, Character s2) {
				// TODO Auto-generated method stub
				return custom.indexOf(s1)-custom.indexOf(s2);
			}
		});
		sum="";
		for(int i=0;i<sum1.length;i++)
			sum += Character.toString(sum1[i]);
		return sum;
	}
	
	public static void add(String r1, String r2){
		String sum = "";
		sum = substitute(r1)+substitute(r2);
		sum = sort(sum);
		System.out.println(reverse_substitute(sum));
		
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String r1,r2;
		Scanner in = new Scanner(System.in);
		System.out.println("Enter first roman numeral");
		r1 = in.nextLine();
		System.out.println("Enter second roman numeral");
		r2 = in.nextLine();
		if(r1.matches("[IVXLCDM]+") && r2.matches("[IXVLCDM]+"))
			add(r1,r2);
		else
			System.out.println("Invalid input. Please enter valid roman letters");
		in.close();
	}
}
