#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author Ben Pritchard
"""
Format the NIST element electron configurations in a table

The NIST web site provides data on the electron configurations for many
elements. This data is not conveniently formatted and this script puts 
that data into a table containing the atom number, electron configuration,
spin multiplicity, and angular momentum, on every line.

To run this script enter:

    format_nistmult.py nistfile

The results are written to standard output.
"""

import sys
import re

if sys.version_info < (3, 6, 0):
   print("This script requires at least Python 3.6.0. You are running:")
   print(sys.version)
   quit(2)

if len(sys.argv) != 2:
   print("Usage: format_nistmult.py nistfile")
   quit(1)


path = sys.argv[1]


atomre = re.compile(r'^(\d+) +(\w+) +([\[\]\w\d]+) +(\d+)(\w)(°?)([\d\/]+)')


# First 5 lines are comments
filelines = [ x.strip() for x in open(path).readlines() ]

for line in filelines[:5]:
  print(line)

curatom = None

for line in filelines[5:]:
    matomre = atomre.match(line)
    if matomre:
      print("{:<5} {:<20} {:<5} {:<5}".format(matomre.group(1), matomre.group(3), matomre.group(4), matomre.group(5)))
