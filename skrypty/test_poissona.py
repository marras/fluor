#!/usr/bin/python

from pylab import *
from numpy import *
from scipy.optimize import leastsq
import scipy
import sys

# To chyba dopasowuje rozkład Poissona do histogramu zliczeń. Nie jestem pewien czy działa.

## Parametric function: 'v' is the parameter vector, 'x' the independent var - funkcja do fitu
#fp = lambda v, x: v[0]*v[1]**x*exp(-v[1])/scipy.misc.factorial(x) # rozklad Poissona, v[0]=C v[1]==lambda
## Normal distribution: C, avg, dev^2
fp = lambda v, x: v[0]*exp(-(x-v[1])**2 / v[2])


## Error function - tak bedzie zawsze! To jest funkcja minimalizowana :)
e = lambda v, x, y: (fp(v,x)-y)

if len(sys.argv)!=2:
	print "Incorrect no. of arguments: give me the histogram file!"
	exit (1)

f = open (sys.argv[1], "rt")

line=f.readline()

x=[]
y=[]

all_lines = f.readlines ()

for line in all_lines:
        l = line.strip().split()
	x.append (int(l[0]))
	y.append (int(l[1]))	#read data from file


## Initial parameter value
v0 = [100.0, 2, 1]
v, success = leastsq(e, v0, args=(x,y), maxfev=10000) ## Fitting

## Plot
def plot_fit():
    print 'Estimated parameters: ', v
#    print 'Real parameters: ', v_real
    X = linspace(min(x),max(x),len(x)*5)
    plot(x,y,'ro', X, fp(v,X))

plot_fit()
show()

