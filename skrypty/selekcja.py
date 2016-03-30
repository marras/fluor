#!/usr/bin/python
from sys import argv

# Skrypt wypisujący z pliku positions.txt kolejne polozenia 1 czasteczki - trzeba najpierw włączyć w .cpp generowanie tego pliku

if len(argv) == 2:
	gr = open ("config.dat")
	for l in gr.readlines():
		a=l.split()
		if len (a) < 3:
			continue
		if a[0] == "NATOMS":
			natoms = int (a[2])

if len(argv) != 2 and len(argv) != 3:
	print "Incorrect number of parameters!\nUsage: %s [which_molecule] -- default config.dat search for natoms" % argv[0]
	print "OR: %s [natoms] [which_molecule]" %argv[0]
	exit(1)

fp = open ("positions.txt")

#natoms = int(argv[1])
which = int(argv[2])

while True:
	for i in range (natoms):
		l = fp.readline().strip().split()
		if int(l[0]) == which:
			print l[1]+" "+l[2]+" "+l[3]
	
