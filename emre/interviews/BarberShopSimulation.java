import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class BarberShopSimulation {
  
  private static class Utils {    
    public static void printStatus(String entity, String status) {
      printStatus(entity, status, false);
    }
    
    public static void printStatus(String entity, String status, boolean error) {
      if(error)
        System.err.println(entity + ": " + status);
      else
        System.out.println(entity + ": " + status);
    }
  }  
  
  private static class Barber implements Runnable {
    private Random random = new Random();
    private BarberShop shop;
    private volatile boolean sleeping = false;
    private static final String ME = "Barber";
    
    public void setBarberShop(BarberShop shop) {
      this.shop = shop;
    }
    
    public void giveHairCut(Customer customer) {
      Utils.printStatus(ME, "Started giving a hair cut to customer " + customer.getId());
      try {
        // randomize hair cut time
        Thread.sleep(random.nextInt(100)); 
      } catch(InterruptedException e) {
        Utils.printStatus(ME, "Whoooaaa I'm interrupted during a hair cut", true);
      }
      Utils.printStatus(ME, "Finished giving a hair cut to customer " + customer.getId());
      customer.doneWitHaircut();
    }
    
    @Override
    public void run() {
      while(shop.isOpen()) {
        Customer c = shop.nextCustomer();
        if(c == null) {
          Utils.printStatus(ME, "Noone is waiting for me. I'll go to sleep");
          sleeping = true;
          synchronized(this) {
            try {
              wait(); 
            } catch(InterruptedException e) {
              Utils.printStatus(ME, "Why am I interrupted?", true);
            }
          }
          sleeping = false;
          Utils.printStatus(ME, "Thanks for waking me up.");    
        } else {
          giveHairCut(c);
        }
      }
      
      Utils.printStatus(ME, "I'm done for the day. Going home...");
    }
    
    public boolean isSleeping() { return sleeping; }
  }
  
  private static class Customer implements Runnable {
    private int id;
    private BarberShop shop;
    private boolean done;
    private Random random = new Random();
    private String me = "Customer";
    
    public Customer(int id, BarberShop shop) {
      this.id = id;
      this.shop = shop;
      this.me += " " + id;
    }

    public int getId() { return id; }
    
    public void doneWitHaircut() {
      Utils.printStatus(me, "Thanks for the haircut, see you next time");
      done = true;
      synchronized(this) {
        notify();
      }
    }
    
    @Override
    public void run() {
      try {
        // randomize customer arrivals to the barber shop
        Thread.sleep(random.nextInt(500));
        Utils.printStatus(me, "I'm going to the barber");
        Thread.sleep(random.nextInt(500));
        if(shop.isOpen()) {
          shop.receive(this); 
        } else {
          Utils.printStatus(me, "Shop is closed. I'll come back later. Or maybe not...");
        }
      } catch(InterruptedException e) {
        Utils.printStatus(me, "Got interrupted before going to the barber shop");
      }  catch(IllegalStateException e) {
        Utils.printStatus(me, "I'm leaving because there are no empty chairs");
        done = true;
      } 
      
      while(!done) {
        synchronized(this) {
          try {
            wait(); 
          } catch(Exception e) {}
        }
      }
    }
    
  }
  
  private static class BarberShop implements Runnable {
    public static final int NUM_CHAIRS = 5;
    public static final int DURATION = 10000;
    private static final String ME = "Shop";
    
    private BlockingQueue<Customer> queue = new ArrayBlockingQueue<Customer>(NUM_CHAIRS);
    private Barber barber;
    private volatile boolean open = false;
    
    public BarberShop() {
      this.barber = new Barber();
      this.barber.setBarberShop(this);
    }
    
    public void receive(Customer customer) throws IllegalStateException {
      if(customer != null) {
        Utils.printStatus(ME, "Customer "  + customer.getId() + " arrived");       
        queue.add(customer);  
        synchronized(barber) {
          if(barber.isSleeping()) {
            Utils.printStatus(ME, "Barber wake up! Somebody needs a hair cut");
            barber.notify(); 
          }
        }
      }
    }
    
    @Override
    public void run() {
      long st = System.currentTimeMillis();
      Utils.printStatus(ME, "We're open now");
      open = true;
      new Thread(barber).start();
      
      while(System.currentTimeMillis() - st < DURATION) {
        try {
          Thread.sleep(100); 
        } catch(InterruptedException e) {}
      }
  
      synchronized(barber) {
        Utils.printStatus(ME, "We're closing...");
        open = false;
        barber.notify();
      }
    }
    
    public Customer nextCustomer() {
      return queue.poll();
    }
    
    public boolean isOpen() { return open; }
  }
  
  public static void main(String... args) throws Exception {
    int numCustomers = 10 + new Random().nextInt(20);
    System.out.println(numCustomers + " customers in this simulation");
    
    BarberShop shop = new BarberShop();
    new Thread(shop).start();
    
    // Using a thread pool for customers because there can be many customers
    ExecutorService pool = Executors.newFixedThreadPool(10);
    for(int i = 1; i <= numCustomers; i++) {
      Customer c = new Customer(i, shop);
      pool.submit(c);
    }   
    pool.shutdown();
    pool.awaitTermination(60, TimeUnit.SECONDS);
    if(pool.isTerminated()) {
      System.out.println("All customers done"); 
    } else {
      System.out.println("Some customers are left hanging");
    }
  }
}