#!/usr/bin/python
''' Implementation of the Z algorithm, based on
http://www.notjustrandom.com/2008/11/01/the-z-algorithm/

May be used in left realigning of indels.
'''

import sys
import optparse

def getZ(s):
    result = {}
 
    l = r = 0
    for k in range(1, len(s)):
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
            zOld = result[kOld]
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
        result[k] = zk
        
    return result

if __name__ == '__main__':
    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()
    
    sequence = args[0]
    print sequence
    print getZ(sequence)
    
    
    