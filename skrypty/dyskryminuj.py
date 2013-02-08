#!/usr/bin/python

import os,sys

if (len(sys.argv) != 2 and len(sys.argv) != 3):
	exit ("Usage: "+sys.argv[0]+" cutoff [filename]")

dyskr = int (sys.argv[1])

if (len(sys.argv)==3): #assume the 3rd param will be the name
	name = sys.argv[2]
else:
	name = "test.flu"

fp = open (name, "r")
f2 = open (name+".dyskr"+`dyskr`, "w")

while 1:
	dane = fp.readline()
	count = int(dane)
	if count > dyskr:
		count = dyskr
	f2.write(str(count)+'\n')
