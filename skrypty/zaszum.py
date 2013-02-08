#!/usr/bin/python

#Adds Poissonian noise to the trajectory (lambda must be given as argument)

import sys, math, time
from random import random as rnd

def rnd_poisson (lamb):
    L =  math.exp(-lamb)
    k = 0
    p = 1
    
    while p > L:
         k = k + 1
         u = rnd() 
	 p = p * u
    
    return k - 1


### Main program ###

if len(sys.argv) != 2: exit("Usage: "+sys.argv[0]+" lambda")

lmb = float(sys.argv[1])
print "Lambda = ",lmb

name = "test.flu"

f2 = open (name+".szum"+"%2.3f" % lmb,"w") #wyrzuc szum do pliku
fp = open (name,"r")
line = fp.readline()
line = "0\n"

while line != "":
	count = int (line.strip())
	count += rnd_poisson (lmb)
#	print `count`
	f2.write(`count`+'\n')
#	time.sleep (0.1)
	line = fp.readline()


		
	

