#!/usr/bin/env python

import sys
import os.path
import struct

def read_crc(fp, size):
  fp.seek(size - 8)
  return struct.unpack('<I', fp.read(4))[0]

with open(sys.argv[1], 'rb') as out:
  out_size = 0

  for arg in sys.argv[2:]:
    size = os.path.getsize(arg)
    out_size += size
    with open(arg, 'rb') as cur:
      cur_crc = read_crc(cur, size)
      out_crc = read_crc(out, out_size)

      m = '%s=%x	%s=%x' % (sys.argv[1], cur_crc, arg, out_crc)
      print m
      if cur_crc != out_crc:
        raise Exception("CRC mistmatch: " + m)
