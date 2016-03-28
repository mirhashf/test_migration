import java.util.*;

/**
 * Created by abhinit on 3/27/16.
 */
public class RomanIntegerAdd {
    /*Below map defines mapping between Roman letters and their value in decimal*/
    private static Map<Character,Integer> symbolMap = new LinkedHashMap();

    /*Below Map defines mapping between  roman number in additive form and corr. subtractive form*/
    private static Map<String,String> additiveToReducedRomanMap = new LinkedHashMap();
    /*Below Map defines mapping between roman number in subtractive  form and corr. additive form*/
    private static Map<String,String> subtractiveToAdditiveRomanMap = new LinkedHashMap();

    static {

        symbolMap.put('M',1000);
        symbolMap.put('D',500);
        symbolMap.put('C',100);
        symbolMap.put('L',50);
        symbolMap.put('X',10);
        symbolMap.put('V',5);
        symbolMap.put('I',1);

        additiveToReducedRomanMap.put("CCCC", "CD");
        additiveToReducedRomanMap.put("DCCCC", "CM");
        additiveToReducedRomanMap.put("XXXX", "XL");
        additiveToReducedRomanMap.put("LXXXX", "XC");
        additiveToReducedRomanMap.put("VIIII", "IX");
        additiveToReducedRomanMap.put("IIII","IV");

        subtractiveToAdditiveRomanMap.put("CD","CCCC");
        subtractiveToAdditiveRomanMap.put("CM","DCCCC");
        subtractiveToAdditiveRomanMap.put("XL","XXXX");
        subtractiveToAdditiveRomanMap.put("XC","LXXXX");
        subtractiveToAdditiveRomanMap.put("IX","VIIII");
        subtractiveToAdditiveRomanMap.put("IV","IIII");

    }

    /**
     * Method replaces substractive Roman notation to addItivie notation.Eg: IV to IIII.
     * @param input
     * @return
     */
    public static StringBuilder replaceSubtractiveNotation(String input){
        StringBuilder result = new StringBuilder(input);
        int index =0;

        for (int i=0;i<result.length();i++){
            for(String key : subtractiveToAdditiveRomanMap.keySet()){
                if(result.indexOf(key,i)>=0) {
                    index = result.indexOf(key,i);
                    result = result.replace(index,index+key.length(),subtractiveToAdditiveRomanMap.get(key));
                    i = index+subtractiveToAdditiveRomanMap.get(key).length();
                    break;
                }
            }
        }

        return result;
    }


    /**
     * Method combines roman number to more compact additive form.
     * @param romanNumerals
     * @return
     */
    private static StringBuilder combineRomanNumeral(StringBuilder romanNumerals){
        Character letter = null;
        int count = 0;
        String result = "";
        int i =romanNumerals.length()-1;

        while (i>=0){
            if(letter==null) {
                count++;
            }
            else if(letter.equals(romanNumerals.charAt(i))){
                count++;
            }
            else if(!letter.equals(romanNumerals.charAt(i)) && count>1){
                result = getRomanNumerals(letter, count);
                i++;
                romanNumerals.replace(i,i+count,result);
                count=1;
            }

            letter = romanNumerals.charAt(i);
            i--;
        }

        if(count>0){
            result = getRomanNumerals(letter,count);
            romanNumerals.replace(0,count,result);
        }

        return romanNumerals;
    }


    /**
     * Method takes romanNumeral and its count and returns Roman number in reduced form.
     * @param romanNumeral
     * @param count
     * @return
     */
    private static String getRomanNumerals(char romanNumeral,int count){
        int sum = symbolMap.get(romanNumeral)*count;
        String result = "";

        for(char letter : symbolMap.keySet()){
            if(sum==0)
                return result;

            while (symbolMap.get(letter)<=sum){
                result = result+letter;
                sum= sum - symbolMap.get(letter);
            }
        }

        return result;
    }


    private static Character[] toCharacterArray(StringBuilder input){
        Character[] result =new Character[input.length()];

        for(int i=0;i<input.length();i++){
            result[i] = input.charAt(i);
        }

        return result;
    }

    /**
     * Method replaces additive Roman number terms to subtractive
     * @param input
     * @return
     */
    private static StringBuilder getReducedForm(StringBuilder input){
        int index = 0;

        for (int i=0;i<input.length();i++){
            for(String key : additiveToReducedRomanMap.keySet()){
                if(input.indexOf(key,i)>=0) {
                    index = input.indexOf(key, i);
                    input = input.replace(index,index+key.length(),additiveToReducedRomanMap.get(key));
                    i = index + additiveToReducedRomanMap.get(key).length();
                    break;
                }
            }
        }

        return input;
    }


    private static class SortRomanNumeralsComparator implements Comparator<Character>{
        @Override
        public int compare(Character lA, Character lB) {
            return symbolMap.get(lB)-symbolMap.get(lA);
        }
    }

    /**
     * Below class defines invalid input exception.
     */
    public static class InvalidInputException extends Exception{
        private static final String errorMessage = "Error: Invalid input";

        public InvalidInputException(){
            super(errorMessage);
        }
    }

    /**
     * Method validates input for valid roman letters and valid subtractive forms.
     * @param num
     * @return
     */
    private static boolean validateInput(String num){
        char prev = '\0';
        String subStr = "";

        for(int i=num.length()-1;i>=0;i--){
            if(symbolMap.get(num.charAt(i))==null)
                return false;

            if(i<num.length()-1 && symbolMap.get(num.charAt(i))<symbolMap.get(prev)) {
                subStr =  num.substring(i,i+2);
                if(subtractiveToAdditiveRomanMap.get(subStr)==null)
                    return false;
            }
            prev = num.charAt(i);
        }

        return true;
    }

    /**
     * Method takes input as two Roman numbers and returns the roman number
     * which is sum of two input roman numbers.
     * @param numA
     * @param numB
     * @return
     */
    public static String addRomans(String numA,String numB) throws InvalidInputException{

        //Validate input
        if(!validateInput(numA) || (!validateInput(numB)))
            throw new InvalidInputException();


        StringBuilder result = replaceSubtractiveNotation(numA)
                .append(replaceSubtractiveNotation(numB));
        Character[] characters = toCharacterArray(result);
        Arrays.sort(characters
                , new SortRomanNumeralsComparator());

        result = new StringBuilder();
        for(Character letter : characters)
                result.append(letter);

        //combine and reduce
        return getReducedForm(combineRomanNumeral(result)).toString();
    }

    public static void main(String[] args){
        String numA = "";
        String numB = "";

        try {
            Scanner scanner = new Scanner(System.in);
            System.out.println("Enter first number");
            numA = scanner.next();
            System.out.println("Enter second number");
            numB = scanner.next();

            String result = addRomans(numA, numB);
            System.out.println("Sum: " + result);
        }
        catch(InvalidInputException ex){
            System.out.println(ex.getMessage());
        }
        catch (Exception ex){
            ex.printStackTrace();
        }

    }
}
