import java.util.HashMap;
import java.util.Map;

public class RomanAddition {

	public static void main(String[] args) {
		if (args.length != 2) {
			System.out.println("Usage: java RomanAddition [RomanNumeral1] [RomanNumeral2]");
			return;
		}
		if (args[0].length() == 0) {
			System.out.println(args[1]);
			return;
		}
		if (args[1].length() == 0) {
			System.out.println(args[0]);
			return;
		}
		String sum = compact(mergeStr(uncompact(args[0]), uncompact(args[1])));
		System.out.println(sum);
	}
	
    // Substitute subtractives in roman numerals
	private static String uncompact(String str) {
		return str.replaceAll("IV", "IIII").replaceAll("IX", "VIIII")
				.replaceAll("XL", "XXXX").replaceAll("XC", "LXXXX")
				.replaceAll("CD", "CCCC").replaceAll("CM", "DCCCC");
	}
	
    // Merge and sort two numerals
	private static String mergeStr(String str1, String str2) {
		Map<Character, Integer> table = new HashMap<Character, Integer>();
		table.put('I', 1);
		table.put('V', 2);
		table.put('X', 3);
		table.put('L', 4);
		table.put('C', 5);
		table.put('D', 6);
		table.put('M', 7);
		
		StringBuilder rs = new StringBuilder();
		int p1 = 0, p2 = 0;
		while (p1 < str1.length() || p2 < str2.length()) {
			int n1 = p1 < str1.length() ? table.get(str1.charAt(p1)) : 0;
			int n2 = p2 < str2.length() ? table.get(str2.charAt(p2)) : 0;
			if (n1 == n2) {
				rs.append(str1.charAt(p1 ++));
				rs.append(str2.charAt(p2 ++));
			} else if (n1 < n2) {
				rs.append(str2.charAt(p2 ++));
			} else {
				rs.append(str1.charAt(p1 ++));
			}
		}
		return rs.toString();
	}

    // Reformat string into Roman numeral 
	private static String compact(String str) {
		str = str.replaceAll("IIIII", "V").replaceAll("VV", "X")
				.replaceAll("XXXXX", "L").replaceAll("LL", "C")
				.replaceAll("CCCCC", "D").replaceAll("DD", "M");
		return str.replaceAll("IIII", "IV").replaceAll("VIIII", "IX")
				.replaceAll("XXXX", "XL").replaceAll("LXXXX", "XC")
				.replaceAll("CCCC", "CD").replaceAll("DCCCC", "CM");
	}
}
