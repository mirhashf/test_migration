package com.bina;

/**
 * Created by muj2 on 2/28/16.
 */
public class Cat {
    public enum Sex{
        MALE, FEMALE
    }

    final String name;
    public Sex sex;
    public double weight;

    public Cat(String name, Sex sex, double weight) {
        this.name = name;
        this.sex = sex;
        this.weight = weight;
    }

    public String getName() {
        return name;
    }

    public Sex getSex() {
        return sex;
    }

    public void setSex(Sex sex) {
        this.sex = sex;
    }

    public double getWeight() {
        return weight;
    }

    public void setWeight(double weight) {
        this.weight = weight;
    }

    @Override
    public String toString() {
        return "Cat{" +
                "name='" + name + '\'' +
                ", sex=" + sex +
                ", weight=" + weight +
                '}';
    }
}
