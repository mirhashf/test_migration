
import java.util.Random;

public class Main {

    public static void main(String[] args) {

        if (args.length == 0) {
            Random randomGenerator = new Random();
            //Happy cases
            for (int i = 0; i < 100; i++) {
                int i1 = randomGenerator.nextInt(10000);
                int i2 = randomGenerator.nextInt(10000);
                String roman1 = RomanNumberUtil.convertRoman(i1);
                String roman2 = RomanNumberUtil.convertRoman(i2);
                String expectedResult = RomanNumberUtil.convertRoman(i1 + i2);
                String actualResult = RomanAdditionUtil.addRomans(roman1, roman2);

                if (expectedResult.equals(actualResult)) {
                    System.out.println("Passed !!! Actual: " + actualResult + " Expected: " + expectedResult);
                } else {
                    System.out.println("Failed !!! Actual: " + actualResult + " Expected: " + expectedResult);
                }
            }
            //UnHappy Cases
            System.out.println("Actual: " + RomanAdditionUtil.addRomans(" ", " ") + " Expected: " + "null");
            System.out.println("Actual: " + RomanAdditionUtil.addRomans("X", " ") + " Expected: " + "null");
        } else if (args.length == 2) {
            System.out.println(RomanAdditionUtil.addRomans(args[0], args[1]));
        }
    }
}

