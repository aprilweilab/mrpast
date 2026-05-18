#!/bin/bash
set -ev

FSC=../../../revisions/fsc28_linux64/fsc28

# Based on testing, the MLE seems to converge before 50 iterations
# (and restarts multiple times before that)
ITER=50
JOBS=20

for rep in $(seq 0 49); do
    mkdir -p rep${rep}
    cd rep${rep}/

    # Cleanup after last time
    rm -f *.obs

    # fsc2 overwrites the .par file, so we start fresh every time.
    cp ../1t12.par.saved 1t12.par 
    cp ../1t12.tpl .
    cp ../1t12.est .

    # Simulate the data from the par file
    ../${FSC} -i 1t12.par -n1 -d -s0 -k 7000000 -r ${rep}
    cp 1t12/*.obs .
    rm -rf 1t12

    # Infer the results
    /usr/bin/time -v ../${FSC} -r ${rep} -t 1t12.tpl -n 100000 -d -e 1t12.est -M -L ${ITER} -q -y 5 -c ${JOBS} 2>&1 | tee rep${rep}.log

    cd ..
done
