package com.gaurav.advanced_validator;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class AdvancedPasswordValidatorTest {

	AdvancedPasswordValidator passwordValidator = new AdvancedPasswordValidator();
	
	@Test
	public void testRule1() {
		String password = "abcde";
		String[] ruleKeys = {"1", "2"};
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "abcdef";
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "abcdeF";
		assertTrue(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "FabcdeF";
		assertTrue(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "FFFFaed";
		assertTrue(passwordValidator.parsePassword(password, ruleKeys));
	}
	
	@Test
	public void testRule2() {
		String password = "abcde";
		String[] ruleKeys = {"1", "3"};
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "abcdef";
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "a2cdeF";
		assertTrue(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "2abcde";
		assertTrue(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "a2e2d4";
		assertTrue(passwordValidator.parsePassword(password, ruleKeys));
	}
	
	@Test
	public void testRule3() {
		String password = "abcde3";
		String[] ruleKeys = {"2", "3"};
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "abcdeD";
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "a2cdeF";
		assertTrue(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "2abcdeRR";
		assertTrue(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "a2e2dDD4";
		assertTrue(passwordValidator.parsePassword(password, ruleKeys));
	}
	
	@Test
	public void testRule4() {
		String password = "abcde";
		String[] ruleKeys = {"1", "2", "3"};
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "";
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
		
		password = null;
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "abcdeD";
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "abcdess2";
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "a2cdeF";
		assertTrue(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "2abcdeRR";
		assertTrue(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "a2e2dDD4";
		assertTrue(passwordValidator.parsePassword(password, ruleKeys));
	}
	
	@Test
	public void testRule5() {
		String password = "aS3gadg";
		String[] ruleKeys = {"2", "8"};
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
	}
	
	@Test
	public void testRule6() {
		String password = "aS3gadgroche";
		String[] ruleKeys = {"1", "2", "9"};
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "roch";
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "rochRe";
		assertTrue(passwordValidator.parsePassword(password, ruleKeys));
		
		password = "Rochd";
		assertFalse(passwordValidator.parsePassword(password, ruleKeys));
	}
}
