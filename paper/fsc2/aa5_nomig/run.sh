#!/bin/bash

# Testing shows the solver converges by 25 iterations
ITER=25
JOBS=20
SEED=42

rm -f *.obs
cp aa5.par.saved aa5.par
../../../revisions/fsc28_linux64/fsc28 -i aa5.par -n1 -d -s0 -k 3000000 -r ${SEED} 2>&1
cp aa5/*.obs .
rm -rf aa5

/usr/bin/time -v ../../../revisions/fsc28_linux64/fsc28 -r ${SEED} -t aa5.tpl -n 100000 -d -e aa5.est -M -L ${ITER} -q -y 5 -c ${JOBS}
