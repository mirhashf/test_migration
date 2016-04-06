// Bowen Wang
// Usage: java RomanAdd arg1 arg2
// arg1 is the first Roman number and arg2 is the second roman number
// The program will generate the sum of arg1 and arg2 then print it to the standard output.
// Example: java RomanAdd CLI LIX
import java.util.*;

class Roman {
	public int[] romanNum;
	public Map<Character, Integer> indexes;
	public Map<Integer, Character> dict;

	public Roman() {
		romanNum = new int[7];
		indexes = generateIndex();
		dict = generateDict();
	}

	public Roman(String num) {
		num = num.toUpperCase();
		romanNum = new int[7];
		indexes = generateIndex();
		dict = generateDict();
		getRoman(num, indexes, romanNum);
	}

	public void getRoman(String num, Map<Character, Integer> indexes, int[] romanNum) {
		char[] charArr = num.toCharArray();
		for(int i=0; i<charArr.length; i++) {
			int index = indexes.get(charArr[i]);
			if(i>0 && index>indexes.get(charArr[i-1])) {
				int lastIndex = indexes.get(charArr[i-1]);
				if(index-lastIndex==1) {
					romanNum[lastIndex] = romanNum[lastIndex]-1+4;
				} else if(index-lastIndex==2) {
					romanNum[lastIndex] = romanNum[lastIndex]-1+4;
					romanNum[lastIndex+1] += 1;
				}
			} else {
				romanNum[index]++;
			}
		}
	}

	public Roman add(Roman a) {
		Roman res = new Roman();
		
		for(int i=0; i<romanNum.length; i++) {
			res.romanNum[i] = romanNum[i]+a.romanNum[i];
		}

		for(int i=0; i<res.romanNum.length-1; i++) {
			if(i%2==0) {
				while(res.romanNum[i]>=5) {
					res.romanNum[i+1]++;
					res.romanNum[i] -= 5;
				}
			} else {
				while(res.romanNum[i]>=2) {
					res.romanNum[i+1]++;
					res.romanNum[i] -= 2;
				}
			}
		}
	
		return res;
	}

	public Map<Character, Integer> generateIndex() {
		Map<Character, Integer> indexes = new HashMap<Character, Integer>();
		indexes.put('I',0);
		indexes.put('V',1);
		indexes.put('X',2);
		indexes.put('L',3);
		indexes.put('C',4);
		indexes.put('D',5);
		indexes.put('M',6);
		return indexes;
	}

	public Map<Integer, Character> generateDict() {
		Map<Integer, Character> dict = new HashMap<Integer, Character>();
		dict.put(0,'I');
		dict.put(1,'V');
		dict.put(2,'X');
		dict.put(3,'L');
		dict.put(4,'C');
		dict.put(5,'D');
		dict.put(6,'M');
		return dict;
	}

	// xix = 19 ix = 9 xiv = 14
	public String getRoman() {
		StringBuilder sb = new StringBuilder();
		for(int i=0; i<romanNum.length; i++) {
			int count = romanNum[i];
			if(count==4 && i<romanNum.length-1) {
				if(romanNum[i+1]>=2 || romanNum[i+1]==0 || i==romanNum.length-2) {
					sb.append(dict.get(i+1));
					sb.append(dict.get(i));
				} else {
					sb.append(dict.get(i+2));
					sb.append(dict.get(i));
					i++;
				}
			} else {
				while(count>0) {
					sb.append(dict.get(i));
					count--;
				}
			}
		}
		String res = sb.reverse().toString();
		return res;
	}

	public void printRoman() {
		System.out.println(getRoman());
	}
}

public class RomanAdd {
	
	public static void main(String[] args) {
		if(args.length!=2) {
			throw new IllegalArgumentException("Usage: Roman1 Roman2");
		}
		if((!checkInput(args[0])) || (!checkInput(args[1]))) {
			throw new IllegalArgumentException("Input can only caontain Roman numeric symbol");
		}

		Roman one = new Roman(args[0]);
		Roman two = new Roman(args[1]);
		Roman res = one.add(two);
		System.out.println(one.getRoman()+" + "+two.getRoman()+" = "+res.getRoman());
	}

	public static boolean checkInput(String num) {
		Set<Character> set = new HashSet<Character>();
		set.add('I');
		set.add('V');
		set.add('X');
		set.add('L');
		set.add('C');
		set.add('D');
		set.add('M');
		for(char c : num.toCharArray()) {
			if(!set.contains(c)) return false;
		}
		return true;
	}	
}