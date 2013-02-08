#!/usr/bin/python

import sys

if len (sys.argv) < 3:
	print "Usage: "+sys.argv[0]+" traj_file threshold"
	exit (1) 

f = open (sys.argv[1], "rt")
fout = open ("on-off.txt", "wt")

prog = float (sys.argv[2])
czas = 1

k = int (f.readline())
if (k > prog): stan = 1  #on
else: stan = 0            #off

while True:
	l = f.readline()
	if not l: break
	k = int (l)
	print k, stan, czas
	if (k > prog and stan == 0 or k < prog and stan == 1):
		fout.write(`czas`+' '+`stan`+'\n')
		czas = 0
		if stan==1: stan = 0
		elif stan==0: stan = 1
		else: print "error!" 
	czas += 1

if czas != 0:
	fout.write(`czas`+' '+`stan`+'\n')

f.close()
fout.close()
