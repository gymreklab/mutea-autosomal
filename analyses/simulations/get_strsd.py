#!/usr/bin/env python

import math
import sys

try:
    p = float(sys.argv[1])
except:
    sys.stderr.write("Usage: ./get_strsd.py <p>\n")
    sys.exit(1)

strsd = math.sqrt((2-p)*1.0/(p**2))
sys.stdout.write(str(strsd)+"\n")
sys.exit(0)
