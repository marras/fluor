#!/usr/bin/python

import sys

f = open (sys.argv[1], "rt")
fout = open ("binnedTraj.flu", "wt")

binL = 50 

while True:
	k = 0;
	for i in range (binL):
		k += int (f.readline())
	fout.write(`k`+'\n')
