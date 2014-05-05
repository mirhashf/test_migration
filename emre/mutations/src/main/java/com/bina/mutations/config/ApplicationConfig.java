package com.bina.mutations.config;

import org.glassfish.jersey.jackson.JacksonFeature;
import org.glassfish.jersey.server.ResourceConfig;

public class ApplicationConfig extends ResourceConfig {

	public ApplicationConfig() {
		packages("com.bina.mutations").register(JacksonFeature.class);
	}
}