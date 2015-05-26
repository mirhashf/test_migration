package com.bina.interview.hangman.logic;

import java.io.PrintStream;
import java.util.Set;

public class SimpleHangmanLogic implements IHangmanLogic {
  private final static String Winner = "winner!";

  private final static String Loser = "loser!";

  private final static char MaskLetter = '_';

  private final int maxBadGuessCount;

  private final Set<Character> hiddenLetters;

  private int badGuessCount;

  public SimpleHangmanLogic(final int maxBadGuessCount, final Set<Character> hiddenLetters) {
    this.maxBadGuessCount = maxBadGuessCount;
    this.hiddenLetters = hiddenLetters;
  }

  @Override
  public boolean simulateTurn(final char letter) {
    if (hiddenLetters.remove(letter)) { return true; }
    badGuessCount += 1;
    return false;
  }

  private boolean isLoser() {
    return badGuessCount >= maxBadGuessCount;
  }

  private boolean isWinner() {
    return hiddenLetters.isEmpty();
  }

  protected String getMaskedWord(final String word) {
    final StringBuilder sb = new StringBuilder();
    for (final char letter : word.toCharArray()) {
      sb.append(hiddenLetters.contains(letter) ? MaskLetter : letter);
    }
    return sb.toString();
  }

  @Override
  public void render(final String word, final PrintStream ps) {
    if (isWinner() || isLoser()) {
      ps.append(String.format("\n%s\n%s", word, isWinner() ? Winner : Loser));
    } else {
      ps.append(String.format("%d %d %s", badGuessCount, hiddenLetters.size(), getMaskedWord(word)));
    }
  }

  @Override
  public boolean isFinished() {
    return isWinner() || isLoser();
  }
}
