package com.gaurav.rules_interface;

@FunctionalInterface
public interface IRule {
	public boolean validate(String password);
}
