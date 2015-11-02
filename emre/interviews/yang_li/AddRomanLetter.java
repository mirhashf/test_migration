/* package whatever; // don't place package name! */

import java.io.*;
import java.util.*;

class myCode
{
    public static HashMap<String, Integer> order;

    public static HashMap<Integer, String> lookUp;

    public static void main (String[] args) throws java.lang.Exception
    {
        order = new HashMap<String, Integer>();
        order.put("I", 1);
        order.put("V", 2);
        order.put("X", 3);
        order.put("L", 4);
        order.put("C", 5);
        order.put("D", 6);
        order.put("M", 7);

        lookUp = new HashMap<Integer, String>();
        lookUp.put(1, "I");
        lookUp.put(2, "V");
        lookUp.put(3, "X");
        lookUp.put(4, "L");
        lookUp.put(5, "C");
        lookUp.put(6, "D");
        lookUp.put(7, "M");

        System.out.println("Calculating 'CCCLXVIIII' + 'DCCCXXXXV': 369 + 845 = 1214");
        System.out.println(addRomanNumeral("CCCLXVIIII", "DCCCXXXXV"));
    }

    public static String addRomanNumeral(String a, String b) {
        a = uncompact(a);
        b = uncompact(b);
        String s = sortRomanNumeral(a + b);
        s = combineSameGroup(s);
        s = compact(s);
        return s;
    }

    // Convert to subtractive forms. For example, IX --> VIIII, IV --> IIII.
    private static String uncompact(String str) {
        String s1;
        String s2;
        StringBuilder sb = new StringBuilder();
        for (int i=0; i<str.length()-1; i++) {
            s1 = str.substring(i, i+1);
            s2 = str.substring(i+1, i+2);
            if (compareTo(s1, s2) < 0) {
                if (order.get(s2) - order.get(s1) == 1) {
                    int index = order.get(s1);
                    for (int j=0; j <4; j++) { sb.append(lookUp.get(index)); }
                } else if (order.get(s2) - order.get(s1) == 2) {
                    int index1 = order.get(s1);
                    int index2 = order.get(s2)-1;
                    sb.append(lookUp.get(index2));
                    for (int j=0; j<4; j++) { sb.append(lookUp.get(index1)); }
                }
                i += 1;
                if (i == str.length()-2) {
                    sb.append(str.substring(i+1, i+2));
                }
            } else {
                sb.append(s1);
                if (i == str.length()-2) {
                    sb.append(s2);
                }
            }
        }
        return sb.toString();
    }

    // Sort ruman numeral from left to right by puting the largest symbol on the left
    private static String sortRomanNumeral(String s) {
        StringBuilder sorted = new StringBuilder();
        String compareS;
        for (int i=7; i>=1; i--) {
            compareS = lookUp.get(i);
            for (int j=0; j<s.length(); j++) {
                if (s.substring(j, j+1).equals(compareS)) {
                    sorted.append(compareS);
                }
            }
        }
        return sorted.toString();
    }

    // compare two strings, return positive (s1>s2), zero (s1=s2), or negative (s1<s2)
    private static int compareTo(String s1, String s2) {
        return order.get(s1) - order.get(s2);
    }

    // Convert the subtractives to the other form. For example, IIII --> IV, VIIII --> IX, etc.
    private static String compact(String s) {
        StringBuilder compacted = new StringBuilder();
        String ret = "";
        String lastNum = s.substring(s.length()-1);
        String currNum = "";
        int repeat = 1;
        for (int i=s.length()-2; i>=0; i--) {
            if (i != -1) { currNum = s.substring(i, i+1); }
            if (lastNum.equals(currNum)) {
                repeat += 1;
            } else {
                if (repeat == 4) {
                    if (order.get(currNum) - order.get(lastNum) == 1) {
                        compacted.append(lastNum);
                        compacted.append(lookUp.get(order.get(currNum)+1));
                        s = s.substring(0, i) + compacted.toString() + s.substring(i+5);
                    } else {
                        compacted.append(lastNum);
                        compacted.append(lookUp.get(order.get(lastNum)+1));
                        s = s.substring(0, i+1) + compacted.toString() + s.substring(i+5);
                    }
                }
                compacted = new StringBuilder();
                lastNum = currNum;
                repeat = 1;
            }
        }
        return s;
    }

    // Combine roman numerals to "larger" symbol, there are two general forms:
    // 1. VV --> X; 2. IIIII --> V.
    private static String combineSameGroup(String s) {
        StringBuilder combined = new StringBuilder();
        String ret = s;
        String lastNum = s.substring(s.length()-1);
        String currNum = "";
        int repeat = 1;
        for (int i=s.length()-2; i>=-1; i--) {
            if (i != -1) { currNum = s.substring(i, i+1); }
            if (i != -1 && lastNum.equals(currNum)) {
                repeat += 1;
            } else {
                if (order.get(lastNum) % 2 == 1 && repeat >= 5) {
                    combined.append(lookUp.get(order.get(lastNum)+1));
                    for (int j=5; j<repeat; j++) { combined.append(lastNum); }
                    ret = ret.substring(0, i+1) + combined.toString() + ret.substring(i+1+repeat);
                    if (i != -1) { i += 1; }
                } else if (order.get(lastNum) % 2 == 0 && repeat >= 2) {
                    combined.append(lookUp.get(order.get(lastNum)+1));
                    for (int j=2; j<repeat; j++) { combined.append(lastNum); }
                    ret = s.substring(0, i+1) + combined.toString() + ret.substring(i+1+repeat);
                    if (i != -1) { i += 1; }
                } else {
                    for (int j=0; j<repeat; j++) { combined.append(lastNum);}
                    ret = ret.substring(0, i+1) + combined.toString() + ret.substring(i+1+repeat);
                }
                combined = new StringBuilder();
                lastNum = currNum;
                repeat = 1;
            }
        }
        return ret;
    }
}