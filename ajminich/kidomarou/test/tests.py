#!/usr/bin/env python
'''Runs unit tests on the Kidomarou pipeline configurator.

Usage: python tests.py
'''

import os
import sys
import unittest
import xmlrunner

# Import modules under test
module_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(module_dir)

import kidomarou

# Universal test parameters
TEST_NAME = "Data Files"
TEST_FILE_TYPE = "test"
IS_MULTIPLE_FILES = True

class TestFileGroup(unittest.TestCase):
    
    def setUp(self):
        self._file_group = kidomarou.FileGroup(
           TEST_NAME, TEST_FILE_TYPE, IS_MULTIPLE_FILES)
        
    def testName(self):
        self.assertEqual(TEST_NAME, self._file_group.get_name())

    def testFileType(self):
        self.assertEqual(TEST_FILE_TYPE, self._file_group.get_file_type())
        
    def testIsMultipleFiles(self):
        self.assertEqual(IS_MULTIPLE_FILES, self._file_group.get_is_multiple_files())
        
class TestFastqGroup(unittest.TestCase):
    
    def setUp(self):
        self._fastq_group = kidomarou.FastqGroup(TEST_NAME, IS_MULTIPLE_FILES)
        
    def testName(self):
        self.assertEqual(TEST_NAME, self._fastq_group.get_name())
    
    def testFileType(self):
        self.assertEqual(kidomarou.FileGroup.FASTQ, self._fastq_group.get_file_type())
        
    def testIsMultipleFiles(self):
        self.assertEqual(IS_MULTIPLE_FILES, self._fastq_group.get_is_multiple_files())
        
    def testGetId(self):
        self.assertEqual(kidomarou.FileGroup.FASTQ + "_multi", 
                         self._fastq_group.get_id())
        
class TestAlignedReadsGroup(unittest.TestCase):
    
    IS_BAM = True
    DIVISION = "by_lane"
    SORT_ORDER = "none"
    
    def setUp(self):
        self._aligned_group = kidomarou.AlignedReadsGroup(TEST_NAME, 
          self.IS_BAM, IS_MULTIPLE_FILES, self.DIVISION, self.SORT_ORDER)

    def testName(self):
        self.assertEqual(TEST_NAME, self._aligned_group.get_name())

    def testFileType(self):
        self.assertEqual(kidomarou.FileGroup.BAM, 
                         self._aligned_group.get_file_type())
        
    def testIsMultipleFiles(self):
        self.assertEqual(IS_MULTIPLE_FILES, 
                         self._aligned_group.get_is_multiple_files())
        
    def testGetId(self):
        self.assertEqual(kidomarou.FileGroup.BAM + "_" + self.DIVISION + 
                         "_" + self.SORT_ORDER, self._aligned_group.get_id())

class TestVariantCallsGroup(unittest.TestCase):
    
    IS_VCF = True
    DIVISION = "by_region"
    
    def setUp(self):
        self._variants_group = kidomarou.VariantCallsGroup(
          TEST_NAME, self.IS_VCF, IS_MULTIPLE_FILES, self.DIVISION)

    def testName(self):
        self.assertEqual(TEST_NAME, self._variants_group.get_name())

    def testFileType(self):
        self.assertEqual(kidomarou.FileGroup.VCF, 
                         self._variants_group.get_file_type())
        
    def testIsMultipleFiles(self):
        self.assertEqual(IS_MULTIPLE_FILES, 
                         self._variants_group.get_is_multiple_files())
        
    def testGetId(self):
        self.assertEqual(kidomarou.FileGroup.VCF + "_" + self.DIVISION + 
                         "_left_aligned", self._variants_group.get_id())

if __name__ == '__main__':
    
    args = sys.argv
    
    if len(args) < 2:
        print __doc__
        sys.exit(1)
    
    # Extract the output folder from the arguments list, and send the remaining arguments
    # to the unit tester
    output_folder = args.pop(1)
    
    unittest.main(argv=args, testRunner=xmlrunner.XMLTestRunner(output=output_folder, verbose=True))
