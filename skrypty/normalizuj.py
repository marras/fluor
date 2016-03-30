#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Normalizuje wszystkie pliki o nazwach corr*.dat (chyba?)

import os,re,sys

def znormalizuj (plik_in, plik_out):
	stary_stdout = sys.stdout
	sys.stdin = plik_in
	sys.stdout = plik_out

	inp = raw_input() #pierwsza linijka to tytul
	print inp
	inp = raw_input()

	inp.strip()
	lista = inp.split()	#split the two numbers

	G0 = float(lista[1])

	try:
		while inp!='':
			inp.strip()
			lista = inp.split()
			lista[1] = float (lista[1])
			lista[1] /= G0
			print lista[0], lista[1]
			inp=raw_input()	

	except EOFError:
		sys.stdout = stary_stdout
		print "Koniec pliku."
		return


dirList = os.listdir (os.getcwd())
for d in dirList:
	if re.match (r'^corr.*\.dat$',d):
		print d, "... ",
		file = open(d,'r')
		out = open (d+".norm",'w')
		znormalizuj (file, out)



