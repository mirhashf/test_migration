package com.gaurav.advanced_validator;

import com.gaurav.password_validator.PasswordValidator;
import com.gaurav.rules_interface.IRule;
import com.gaurav.rules_validator_interface.IPasswordValidator;

public class AdvancedPasswordValidator implements IPasswordValidator{
	PasswordValidator passwordValidator;
	public AdvancedPasswordValidator() {
		passwordValidator = new PasswordValidator();
		defineRules();
	}

	public void defineRules() {
		passwordValidator.addRule(9, (password) -> {
			if(password == null) return false;
			return !(password.indexOf("bina") != -1 || password.indexOf("roche") != -1);
		});
	}

	@Override
	public boolean addRule(Integer ruleNumber, IRule rule) {
		return passwordValidator.addRule(ruleNumber, rule);
	}

	@Override
	public boolean parsePassword(String password, String[] keys) {
		return passwordValidator.parsePassword(password, keys);
	}

	@Override
	public void deleteRule(Integer ruleNumber) {
		passwordValidator.deleteRule(ruleNumber);
	}
	
	public static void main(String[] args) {
		if(args.length != 2) {
			System.out.println("Error: Proper usage is - java AdvancedPasswordValidator [password] [rule 1, rule2 ..]");
			System.exit(1);
		}
		String password = args[0];
		String[] ruleKeys = args[1].split(",");
		AdvancedPasswordValidator validator = new AdvancedPasswordValidator();
		System.out.println(validator.parsePassword(password, ruleKeys));
	}
}
