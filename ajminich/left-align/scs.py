#!/usr/bin/python
''' Determines the shortest completely repetitive sequence in the
string. Used for left aligning insertions and deletions.

Uses an implementation of the Z algorithm based on
http://www.notjustrandom.com/2008/11/01/the-z-algorithm/
'''

import sys
import optparse

def getShortestCommonSequence(s):
    ''' Gets the shortest sequence that is completely repetitive
    across the entire string. Several examples:
    
    actactact => shortest common sequence is act.
    actgtca   => shortest common sequence is actgtca.
    aaaaaaaaa => shortest common sequence is a.
    
    Utilizes a modified Z algorithm.
    '''
    
    result = []
 
    # The longest possible common sequence is the entire string.
    # The second longest possible sequence is half of the string,
    # so we will only check up to half of the string before
    # simply returning the entire string.
    maxCommonSeqLength = len(s) / 2 + 1
 
    l = r = 0
    for k in range(1, maxCommonSeqLength):
        if k > r:
            zk = 0
            for si in range(0, len(s)):
                if k + si < len(s) and \
                    s[si] == s[k + si]:
                    pass
                else:
                    break
            if si > 0:
                zk = si
                r = zk + k - 1
                l = k
        else:
            kOld = k - l
            zOld = result[kOld - 1]
            b = r - k + 1
            if zOld < b:
                zk = zOld
            else:
                zk = b
 
                for si in range(b, len(s)):
                    if k + si < len(s) and s[si] == s[k + si]:
                        pass
                    else:
                        break
                zk = si
                r = zk + k - 1
                l = k
        result.append(zk)
        
        # Check if the length of P is actually the rest of the string
        if zk == len(s) - k:
            return s[0:k]
        
    return s

if __name__ == '__main__':
    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()
    
    sequence = args[0]
    
    scs = getShortestCommonSequence(sequence)
    
    print "Input sequence is '" + sequence + "'."
    print "Shortest common sequence is '" + scs + "'."
    
    
    