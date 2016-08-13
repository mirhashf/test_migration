package com.gaurav.password_validator;

import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;

import com.gaurav.rules_interface.IRule;
import com.gaurav.rules_validator_interface.IPasswordValidator;

public final class PasswordValidator implements IPasswordValidator {
	private static class PasswordRule {
		private final static Map<Integer, IRule> password_rules = new HashMap<>();
	}

	public PasswordValidator() {
		defineRules();
	}

	public void defineRules() {
		PasswordRule.password_rules.put(1, (password) -> {
			if(password == null) return false;
			return password.length() >= 6;
		});

		PasswordRule.password_rules.put(2, (password) -> {
			if(password == null) return false;
			final Pattern oneUpperCase = Pattern.compile("^(?=.*[A-Z]).+$");
			return oneUpperCase.matcher(password).matches();
		});

		PasswordRule.password_rules.put(3, (password) -> {
			if(password == null) return false;
			final Pattern oneNumericChar = Pattern.compile("^(?=.*\\d).+$");
			return oneNumericChar.matcher(password).matches();
		});
	}

	@Override
	public boolean addRule(Integer ruleNumber, IRule rule) {
		boolean ruleAdded = false;
		if(!PasswordRule.password_rules.containsKey(ruleNumber)) {
			PasswordRule.password_rules.put(ruleNumber, rule);
			ruleAdded = true;
		}
		return ruleAdded;
	}

	@Override
	public boolean parsePassword(String password, String[] rules) {
		try {
			for(String rule: rules) {
				if(!PasswordRule.password_rules.containsKey(Integer.valueOf(rule)) || 
						!PasswordRule.password_rules.get(Integer.valueOf(rule)).validate(password)) {
					return false;
				}
			}
		} catch (NumberFormatException e) {
			System.out.println("Error: Please enter valid rule numbers");
			System.exit(1);
		}
		return true;
	}

	@Override
	public void deleteRule(Integer ruleNumber) {
		if(PasswordRule.password_rules.containsKey(ruleNumber))
			PasswordRule.password_rules.remove(ruleNumber);
	}

	public static void main(String[] args) {
		if(args.length != 2) {
			System.out.println("Error: Proper usage is - java PasswordValidator [password] [rule 1, rule2 ..]");
			System.exit(1);
		}
		String password = args[0];
		String[] ruleKeys = args[1].split(",");
		PasswordValidator passwordValidator = new PasswordValidator();
		System.out.println(passwordValidator.parsePassword(password, ruleKeys));
	}
}
