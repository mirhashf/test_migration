package com.gaurav.rules_validator_interface;

import com.gaurav.rules_interface.IRule;

public interface IPasswordValidator {
	public boolean parsePassword(String password, String[] ruleNumber);
	public boolean addRule(Integer ruleNumber, IRule rule);
	public void deleteRule(Integer ruleNumber);
}
