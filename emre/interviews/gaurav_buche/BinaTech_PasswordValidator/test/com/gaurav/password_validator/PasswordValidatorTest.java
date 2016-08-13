package com.gaurav.password_validator;

import org.junit.Test;
import static org.junit.Assert.*;


public class PasswordValidatorTest {
	PasswordValidator passwordValidator = new PasswordValidator();
	
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
}
