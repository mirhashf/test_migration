// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.util.HashMap;
import java.util.Map;
import java.util.Properties;

import org.springframework.beans.BeansException;
import org.springframework.beans.factory.config.ConfigurableListableBeanFactory;
import org.springframework.beans.factory.config.PropertyPlaceholderConfigurer;


/**
 * Reads properties from the spring.properties file.
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 */
public class SpringPropertiesUtil extends PropertyPlaceholderConfigurer {

  private static Map<String, String> propertiesMap;
  // Default as in PropertyPlaceholderConfigurer
  private int springSystemPropertiesMode = SYSTEM_PROPERTIES_MODE_FALLBACK;

  @Override
  public void setSystemPropertiesMode(int systemPropertiesMode) {
    super.setSystemPropertiesMode(systemPropertiesMode);
    springSystemPropertiesMode = systemPropertiesMode;
  }

  @Override
  protected void processProperties(ConfigurableListableBeanFactory beanFactory, Properties props)
      throws BeansException {
    super.processProperties(beanFactory, props);
    propertiesMap = new HashMap<String, String>();
    for (Object key : props.keySet()) {
      String keyStr = key.toString();
      String valueStr = resolvePlaceholder(keyStr, props, springSystemPropertiesMode);
      propertiesMap.put(keyStr, valueStr);
    }
  }

  public static String getProperty(String name) {
    return propertiesMap.get(name).toString();
  }

}
