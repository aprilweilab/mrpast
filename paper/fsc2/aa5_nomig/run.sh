#!/bin/bash

# Testing shows the solver converges by 25 iterations
ITER=25
JOBS=20
SEED=42

DATASET=../fsc2_aa5_20full

rm -f *.obs
cp ${DATASET}/*.obs .
echo "Using SFS from ${DATASET}"

FSC=../../../../revisions/fsc28_linux64/fsc28

# Cleanup old files and copy the model here.
rm -f aa5.par
cp ../aa5.tpl .
cp ../aa5.est .

# We use "-0" because I don't know how robust the non-segregating site info from msprime is,
# it likely does not help/match the assumptions of fsc2
/usr/bin/time -v ${FSC} -r ${SEED} -t aa5.tpl -0 -n 100000 -d -e aa5.est -M -L ${ITER} -q -y 5 -c ${JOBS}
