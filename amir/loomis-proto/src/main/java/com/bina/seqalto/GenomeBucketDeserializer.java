// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.io.IOException;

import org.codehaus.jackson.JsonProcessingException;
import org.codehaus.jackson.map.DeserializationContext;
import org.codehaus.jackson.map.KeyDeserializer;
import org.codehaus.jackson.map.ObjectMapper;

/**
 * Deserializer to convert {@link GenomeBucket} to JSON.
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 * Created on Feb 13, 2012
 */
public class GenomeBucketDeserializer extends KeyDeserializer{
  /* (non-Javadoc)
   * @see org.codehaus.jackson.map.KeyDeserializer#deserializeKey(java.lang.String, org.codehaus.jackson.map.DeserializationContext)
   */
  @Override
  public Object deserializeKey(String key, DeserializationContext ctxt) throws IOException,
      JsonProcessingException {
    
    ObjectMapper mapper = new ObjectMapper();
    
    return mapper.readValue(key,GenomeBucket.class);
  }
}
