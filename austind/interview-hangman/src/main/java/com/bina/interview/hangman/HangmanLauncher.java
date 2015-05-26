package com.bina.interview.hangman;

import java.io.IOException;

import com.bina.interview.hangman.hidden.IHiddenLetterGenerator;
import com.bina.interview.hangman.hidden.ShuffleLetterGenerator;

public class HangmanLauncher {
  /*
   * java -jar target/interview-hangman.jar <max guess count> <missing letter count> <word>
   */
  public static void main(final String[] args) throws IOException {
    if (args.length != 3) { throw new IllegalArgumentException(); }
    HangmanLauncher.create(args[0], args[1], args[2], new ShuffleLetterGenerator()).mainLoop();
  }

  protected static HangmanRunner create(final String _maxGuessCount, final String _missingLetterCount, final String word, final IHiddenLetterGenerator generator) {
    final int maxBadGuessCount = Integer.parseInt(_maxGuessCount);
    if (maxBadGuessCount < 0) { throw new IllegalArgumentException(); }
    final int missingLetterCount = Integer.parseInt(_missingLetterCount);
    if (missingLetterCount <= 0) { throw new IllegalArgumentException(); }
    if (word.contains("_")) { throw new IllegalArgumentException(); }
    if (StringUtil.toSet(word).size() < missingLetterCount) { throw new IllegalArgumentException(); }
    return new HangmanRunner(maxBadGuessCount, word, generator.generate(word, missingLetterCount));
  }
}
