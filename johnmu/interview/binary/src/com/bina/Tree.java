package com.bina;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by muj2 on 2/28/16.
 */
public class Tree {
    protected List<Cat> list = new ArrayList<>();

    public Tree(List<Cat> list) {
        this.list = list;
    }

    public Tree() {
    }

    public void addCat(Cat c) {
        list.add(c);
    }

    public Cat findByWeight(final double weight) {
        sortListByWeight();
        int lowIndex = 0;
        int highIndex = list.size() - 1;
        int middle = (int) Math.floor((highIndex - lowIndex) / 2);
        while (highIndex != lowIndex) {
            if (list.get(middle).getWeight() > weight) {
                highIndex = middle - 1;
                middle = (int) Math.floor((highIndex - lowIndex) / 2);
            } else {
                lowIndex = middle + 1;
                middle = (int) Math.floor((highIndex - lowIndex) / 2);
            }
        }
        return list.get(highIndex);
    }

    public void sortListByWeight() {
        // Bubble sort FTW
        boolean swapped = true;
        while (swapped) {
            swapped = false;
            for (int i = 0; i < list.size() - 1; i++) {
                if (list.get(i).getWeight() > list.get(i + 1).getWeight()) {
                    Cat temp = list.get(i);
                    list.set(i, list.get(i + 1));
                    list.set(i + 1, temp);
                    swapped = true;
                }
            }
        }
    }

    public double averageWeight(Cat.Sex sex) {
        double total = 0;
        int count = 0;
        for (Cat c : list) {
            if (c.sex == sex) {
                total += c.getWeight();
                count++;
            }
        }
        return total/count;
    }

    public Collection<String> getNames(){
        return list.stream().map((c)->c.name).collect(Collectors.toList());
    }
}
