package com.bina.interview.hangman;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Set;

import com.bina.interview.hangman.logic.IHangmanLogic;
import com.bina.interview.hangman.logic.SimpleHangmanLogic;

public class HangmanRunner {
  private static final String Smile = ":)";

  private static final String Frown = ":(";

  private final String word;

  private final IHangmanLogic state;

  public HangmanRunner(final int maxBadGuessCount, final String word, final Set<Character> hiddenLetters) {
    state = new SimpleHangmanLogic(maxBadGuessCount, hiddenLetters);
    this.word = word;
  }

  public void mainLoop() throws IOException {
    try (final BufferedReader reader = new BufferedReader(new InputStreamReader(System.in))) {
      while (true) {
        state.render(word, System.out);
        System.out.println();
        if (state.isFinished()) {
          break;
        }
        System.out.print("> ");
        System.out.println();
        final String line = reader.readLine();
        if (line == null) {
          break;
        } else if (line.length() != 1) {
          continue;
        }
        System.out.print(state.simulateTurn(line.charAt(0)) ? Smile : Frown);
        System.out.println();
      }
    }
  }
}
