#! /usr/bin/python

# A Python implementation of the forward and inverse Burrows-Wheeler Transform.
# taken from http://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform.

import sys

def bwt(s):
    """Apply Burrows-Wheeler transform to input string."""
    assert "\0" not in s, "Input string cannot contain null character ('\0')"
    s += "\0"  # Add end of file marker
    table = sorted(s[i:] + s[:i] for i in range(len(s)))  # Table of rotations of string
    last_column = [row[-1:] for row in table]  # Last characters of each row
    return "".join(last_column)  # Convert list of characters into string
    
def ibwt(r):
    """Apply inverse Burrows-Wheeler transform."""
    table = [""] * len(r)  # Make empty table
    for i in range(len(r)):
        table = sorted(r[i] + table[i] for i in range(len(r)))  # Add a column of r
    s = [row for row in table if row.endswith("\0")][0]  # Find the correct row (ending in "\0")
    return s.rstrip("\0")  # Get rid of trailing null character
    
if __name__ == '__main__':
    
    input = sys.argv[1]
    
    print "Input is: %s" % input
    
    t = bwt(input)
    
    print "BWT is: %s" % t
    
    it = ibwt(t)
    
    print "IBWT is: %s" % it