#!/bin/bash

A=`pwd | sed -e 's/^.*eps_//'`

echo $A

#./fluor -f
./fluor -b
mv correlationLOG.dat "corr_$A.dat"
./dyskryminuj.py 1
mv test.flu.dyskr test.flu
./fluor -b
mv correlationLOG.dat "corr_${A}_Dyskr1.dat"
