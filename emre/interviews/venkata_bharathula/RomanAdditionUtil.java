
/**
 * Created by venkata on 10/20/15.
 */


import java.util.*;

public final class RomanAdditionUtil {


    private RomanAdditionUtil() {
    }

    private static String compressRoman(String roman) {
        Map<String, String> strMap = new HashMap<String, String>();
        strMap.put("IIIII", "V");
        strMap.put("VV", "X");
        strMap.put("XXXXX", "L");
        strMap.put("LL", "C");
        strMap.put("CCCCC", "D");
        strMap.put("DD", "M");
        for (Map.Entry<String, String> pair : strMap.entrySet()) {
            String key;
            key = pair.getKey();
            String value;
            value = pair.getValue();
            roman = roman.replaceAll(key, value);
        }
        return roman;
    }


    private static String remapRoman(String roman) {
        Map<String, String> strMap = new LinkedHashMap<String, String>();
        strMap.put("DCCCC", "CM");
        strMap.put("CCCC", "CD");
        strMap.put("LXXXX", "XC");
        strMap.put("XXXX", "XL");
        strMap.put("VIIII", "IX");
        strMap.put("IIII", "IV");

        for (Map.Entry<String, String> pair : strMap.entrySet()) {
            String key;
            key = pair.getKey();
            String value;
            value = pair.getValue();
            roman = roman.replaceAll(key, value);
        }
        return roman;
    }

    private static String mergeRoman(String roman1, String roman2) {
        HashMap<Character, Integer> strMap = new HashMap<Character, Integer>();
        strMap.put('I', 0);
        strMap.put('V', 1);
        strMap.put('X', 2);
        strMap.put('L', 3);
        strMap.put('C', 4);
        strMap.put('D', 5);
        strMap.put('M', 6);
        StringBuilder db = new StringBuilder();
        int i = 0, j = 0;
        while (i < roman1.length() && j < roman2.length()) {
            if (strMap.get(roman1.charAt(i)) >= strMap.get(roman2.charAt(j))) {
                db.append(roman1.charAt(i));
                i++;
            } else {
                db.append(roman2.charAt(j));
                j++;
            }
        }
        if (j < roman2.length()) db.append(roman2.substring(j, roman2.length()));
        if (i < roman1.length()) db.append(roman1.substring(i, roman1.length()));
        return db.toString();
    }

    private static String expandRoman(String roman) {
        HashMap<String, String> strMap = new HashMap<String, String>();
        strMap.put("IV", "IIII");
        strMap.put("IX", "VIIII");
        strMap.put("XL", "XXXX");
        strMap.put("XC", "LXXXX");
        strMap.put("CD", "CCCC");
        strMap.put("CM", "DCCCC");
        StringBuilder db = new StringBuilder();
        int i;
        for (i = 0; i < roman.length(); i++) {
            if (i < roman.length() - 1 && strMap.containsKey(roman.substring(i, i + 2))) {
                db.append(strMap.get(roman.substring(i, i + 2)));
                i++;
            } else {
                db.append(roman.charAt(i));
            }
        }
        return db.toString();
    }

    private static boolean isRoman(String roman) {
        Set<Character> romanSet = new HashSet<Character>();
        romanSet.add('I');
        romanSet.add('V');
        romanSet.add('X');
        romanSet.add('L');
        romanSet.add('C');
        romanSet.add('D');
        romanSet.add('M');

        for (int i=0; i < roman.length(); i++) {
            if (!romanSet.contains(roman.charAt(i))) return false;
        }
        return true;
    }

    public static String addRomans(String roman1, String roman2) {

        if (!(isRoman(roman1) && isRoman(roman2))) return null;
        String s = mergeRoman(expandRoman(roman1), expandRoman(roman2));
        while (!s.equals(compressRoman(s))) {
            s = compressRoman(s);
        }
        while (!s.equals(remapRoman(s))) {
            s = remapRoman(s);
        }
        return s;
    }

}
