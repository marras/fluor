#!/bin/sh

make
for katalog in sym*; do
 cp fluor $katalog
# cp skrypt_rods.sh $katalog
 cp grompp.mdp $katalog
done

