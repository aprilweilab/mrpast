#!/bin/bash
set -ev

FSC=../../../../../revisions/fsc28_linux64/fsc28

# Based on testing, the MLE seems to converge before 50 iterations
# (and restarts multiple times before that)
ITER=50
JOBS=20

for rep in $(seq 0 49); do
    mkdir -p rep${rep}
    cd rep${rep}/

    # Cleanup after last time
    rm -f *.obs
    rm -rf 1t12/

    cp ../../1t12.tpl .
    cp ../../1t12.est .

    # We already generated GRGs that have bi-allelic SNPs from all 10 chromosomes
    python ../../../grapp_sfs.py ../../simdata_1t12/${rep}/1t12_merged.biallelic.grg 1t12 fsc2_grapp_data

    cp fsc2_grapp_data/*.obs .

    # Infer the results (use -0 to avoid non-segregating sites, since I'm not sure how robust those are)
    /usr/bin/time -v ${FSC} -r ${rep} -t 1t12.tpl -0 -n 100000 -d -e 1t12.est -M -L ${ITER} -q -y 5 -c ${JOBS} 2>&1 | tee rep${rep}.log

    cd ..
done
