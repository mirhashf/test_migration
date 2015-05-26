package com.bina.interview.hangman.logic;

import java.io.PrintStream;

public interface IHangmanLogic {
  public boolean simulateTurn(final char letter);

  public void render(final String word, final PrintStream ps);

  public boolean isFinished();
}
