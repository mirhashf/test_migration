package com.bina;

public class Main {

    public static void main(String[] args) {
	    Tree tree = new Tree();

        // Add some cats
        tree.addCat(new Cat("bar", Cat.Sex.MALE, 100/3.0));
        tree.addCat(new Cat("baz", Cat.Sex.FEMALE, 30/3.0));
        tree.addCat(new Cat("goo", Cat.Sex.FEMALE, 10));
        tree.addCat(new Cat("foo", Cat.Sex.MALE, 200/3));
        tree.addCat(new Cat("gah", Cat.Sex.MALE, 2/3));
        tree.addCat(new Cat("woo", Cat.Sex.MALE, 0/3));
        System.out.println(tree.findByWeight(10));
        System.out.println(tree.findByWeight(33));
        System.out.println(tree.findByWeight(50));
        System.out.println(tree.averageWeight(Cat.Sex.MALE));
    }
}
