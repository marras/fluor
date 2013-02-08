#!/usr/bin/python

from sys import argv
from scipy import *
from pylab import *

if len(argv) < 2: print argv[0],"file bins"; exit(0)

fin = open (argv[1],"rt")
b = int (argv[2])

a = [float(l.split()[0]) for l in fin.readlines()]

n, bins, patches = hist(a, bins=b)
print n, bins, patches

xlabel('State persistence time')
ylabel('Probability')
title(r'$\mathrm{Histogram\ of\ state\ occupation\ times:} \bar\tau = $ unknown')
#axis([40, 160, 0, 0.03])
grid(True)

show()
