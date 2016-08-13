package a;

public class A {
  public static void main(final String[] args) throws ClassNotFoundException {
    Class.forName("b.B");
    Class.forName("c.C");
  }
}
