package com.bina.interview.hangman.hidden;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import com.bina.interview.hangman.StringUtil;

public class ShuffleLetterGenerator implements IHiddenLetterGenerator {
  @Override
  public Set<Character> generate(final String word, final int missingLetterCount) {
    final List<Character> letters = new ArrayList<>(StringUtil.toSet(word));
    Collections.shuffle(letters, new Random(0));
    final Set<Character> missingLetters = new HashSet<>();
    for (int i = 0; i < missingLetterCount; i++) {
      missingLetters.add(letters.get(i));
    }
    return missingLetters;
  }
}
