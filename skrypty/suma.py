#!/usr/bin/python

import sys

f = open (sys.argv[1], "rt")

print sum ([int(l.split(" ")[0]) for l in f.readlines()])
