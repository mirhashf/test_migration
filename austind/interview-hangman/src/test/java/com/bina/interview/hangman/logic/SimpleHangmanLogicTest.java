package com.bina.interview.hangman.logic;

import java.util.Arrays;
import java.util.HashSet;

import org.junit.Assert;
import org.junit.Test;
public class SimpleHangmanLogicTest {
  @Test
  public void testGetMaskedWord() {
    final SimpleHangmanLogic logic = new SimpleHangmanLogic(0, new HashSet<>(Arrays.asList('a', 'b')));
    Assert.assertEquals("__c__d", logic.getMaskedWord("abcabd"));
  }
}
