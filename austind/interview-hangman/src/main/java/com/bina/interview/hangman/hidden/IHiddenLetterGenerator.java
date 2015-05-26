package com.bina.interview.hangman.hidden;

import java.util.Set;

public interface IHiddenLetterGenerator {
  public Set<Character> generate(final String word, final int hiddenLetterCount);
}
