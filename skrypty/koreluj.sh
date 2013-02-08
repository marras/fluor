#!/bin/bash

R=$RANDOM

for dir in "$@"
do
	cat $dir/corr.dat >> cor_all_$R
#	cat $dir/corr.dat.norm >> cor_all_$R
	echo $dir
done

echo $R

xmgrace cor_all_$R &
sleep 3
rm cor_all_$R
