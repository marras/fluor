#!/usr/bin/python
import os, sys

plik1 = open (sys.argv[1],"r")
plik2 = open (sys.argv[2],"r")

while (1):
	a1 = plik1.readline().split()[1]
	a2 = plik2.readline().split()[1]

	srednia = (a1+a2)/2

	print a1, a2, srednia



print liczby1


