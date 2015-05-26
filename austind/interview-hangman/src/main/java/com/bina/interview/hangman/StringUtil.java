package com.bina.interview.hangman;

import java.util.HashSet;
import java.util.Set;

public class StringUtil {
  public static Set<Character> toSet(final String value) {
    final Set<Character> letters = new HashSet<>();
    for (final char letter : value.toCharArray()) {
      letters.add(letter);
    }
    return letters;
  }
}
