// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Stack;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Reads a FASTA index file and comes up with semi equal genomic regions.
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com) Created on Feb 8, 2012
 */
public class GenomeRegionSplitter {
  private static Logger logger = LoggerFactory.getLogger(GenomeRegionSplitter.class);
  private int numberOfBuckets;
  private File fastaFileIndex;

  /**
   * Constructor for a Splitter
   * 
   * @param fastaFile FASTA File to use to split regions on.
   * @param numberOfBuckets
   * @throws FileNotFoundException
   */
  public GenomeRegionSplitter(File fastaFileIndex, int numberOfBuckets)
      throws FileNotFoundException {
    logger.debug("Created GenomeRegionSplitter for FASTA file: " + fastaFileIndex);
    this.numberOfBuckets = numberOfBuckets;
    this.fastaFileIndex = fastaFileIndex;

  }

  public GenomeBucket[] getRegions() throws NumberFormatException, IOException {
    // go over all the regions adding the sizes
    long totalSize = 0;

    BufferedReader reader = new BufferedReader(new FileReader(fastaFileIndex));
    String line = null;


    Stack<GenomeRegion> regions = new Stack<GenomeRegion>();

    while ((line = reader.readLine()) != null) {
      String[] entry = line.split("\t", 3);
      GenomeRegion region = new GenomeRegion(entry[0], 1, Integer.valueOf(entry[1]));
      regions.add(region);
      totalSize += region.getSize();
    }

    // find optimal bucket size
    long optimalRegionSize =  (long) Math.ceil (Double.valueOf(totalSize) / Double.valueOf(numberOfBuckets));


    GenomeBucket[] splitRegions = new GenomeBucket[numberOfBuckets];

    // output data structure

    for (int i = 0; i < numberOfBuckets; i++) {
      splitRegions[i] = new GenomeBucket();
    }


    int bucketId = 0;
    // keep adding the
    while (!regions.isEmpty()) {
      GenomeRegion region = regions.pop();
      long availableCapacity = optimalRegionSize - splitRegions[bucketId].getSize();
      if (region.getSize() >= availableCapacity) {

        // region too big, need to break.
        GenomeRegion subRegion =
            new GenomeRegion(region.getReference(), region.getStart(), region.getStart()
                + availableCapacity);

        splitRegions[bucketId].addRegion(subRegion);
        region.setStart(region.getStart() + availableCapacity);
        if (region.getSize() != 0){
          regions.push(region);
        }
        bucketId++;
      } else {
        splitRegions[bucketId].addRegion(region);
      }
    }

    return splitRegions;
  }



}
