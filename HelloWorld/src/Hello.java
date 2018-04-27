/**
 * Created by dev on 04/25/18
 */
public class Hello {
    public static void main(String[] args){

        System.out.println("Hello, World!");

        int myFirstNumber = (10 + 5) + (2 * 10);
        int mySecondNumber = 12;
        int myThirdNumber = myFirstNumber * 2;

        int myTotal = myFirstNumber + mySecondNumber + myThirdNumber;
        int myLastOne = 1000 - myTotal;

        System.out.println("Your total value is: " + myTotal + ", and the difference of a 1000 is: " + myLastOne);

    }
}
