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
TEST_FILE_TYPE = "test"
IS_MULTIPLE_FILES = True

class TestFileGroup(unittest.TestCase):
    
    def setUp(self):
        self._file_group = kidomarou.FileGroup(TEST_FILE_TYPE, IS_MULTIPLE_FILES)

    def testFileType(self):
        self.assertEqual(TEST_FILE_TYPE, self._file_group.get_file_type())
        
    def testIsMultipleFiles(self):
        self.assertEqual(IS_MULTIPLE_FILES, self._file_group.get_is_multiple_files())
        
class TestFastqGroup(unittest.TestCase):
    
    def setUp(self):
        self._fastq_group = kidomarou.FastqGroup(IS_MULTIPLE_FILES)
    
    def testFileType(self):
        self.assertEqual(kidomarou.FileGroup.FASTQ, self._fastq_group.get_file_type())
        
    def testIsMultipleFiles(self):
        self.assertEqual(IS_MULTIPLE_FILES, self._fastq_group.get_is_multiple_files())
        
    def testGetName(self):
        self.assertEqual(kidomarou.FileGroup.FASTQ + "_multi", 
                         self._fastq_group.get_name())
        
class TestAlignedReadsGroup(unittest.TestCase):
    
    IS_BAM = True
    DIVISION = "by_lane"
    SORT_ORDER = "none"
    
    def setUp(self):
        self._aligned_group = kidomarou.AlignedReadsGroup(
          self.IS_BAM, IS_MULTIPLE_FILES, self.DIVISION, self.SORT_ORDER)

    def testFileType(self):
        self.assertEqual(kidomarou.FileGroup.BAM, 
                         self._aligned_group.get_file_type())
        
    def testIsMultipleFiles(self):
        self.assertEqual(IS_MULTIPLE_FILES, 
                         self._aligned_group.get_is_multiple_files())
        
    def testGetName(self):
        self.assertEqual(kidomarou.FileGroup.BAM + "_" + self.DIVISION + 
                         "_" + self.SORT_ORDER, self._aligned_group.get_name())

if __name__ == '__main__':
    
    args = sys.argv
    
    if len(args) < 2:
        print __doc__
        sys.exit(1)
    
    # Extract the output folder from the arguments list, and send the remaining arguments
    # to the unit tester
    output_folder = args.pop(1)
    
    unittest.main(argv=args, testRunner=xmlrunner.XMLTestRunner(output=output_folder, verbose=True))
