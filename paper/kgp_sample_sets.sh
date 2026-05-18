#!/bin/bash
set -ev

# Build the lists of sample identifiers, grouped by population, for each of our
# datasets. The order here should match the order in the corresponding .popmap.json

# MXL only has 35 unrelated individuals, which is our minimum count, so we use
# this for every population.
PER_POP=35

# Make the sample sets for each of the three admixture datasets, and the OOA3 dataset as well

# OOA3: YRI_CEU_CHB
{ head -n ${PER_POP} samples_by_pop/samples.yri.nofam.txt && 
  head -n ${PER_POP} samples_by_pop/samples.ceu.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.chb.nofam.txt; } > samples.ooa3.i35.txt

# ADMIX_PUR: YRI_CEU_CHB_PUR
{ head -n ${PER_POP} samples_by_pop/samples.yri.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.ceu.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.chb.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.pur.nofam.txt; } > samples.admix_pur.i35.txt

# ADMIX_CLM: YRI_CEU_CHB_CLM
{ head -n ${PER_POP} samples_by_pop/samples.yri.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.ceu.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.chb.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.clm.nofam.txt; } > samples.admix_clm.i35.txt

# ADMIX_MXL: YRI_CEU_CHB_MXL
{ head -n ${PER_POP} samples_by_pop/samples.yri.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.ceu.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.chb.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.mxl.nofam.txt; } > samples.admix_mxl.i35.txt

# ADMIX_CLM: YRI_GBR_CHB_CLM
{ head -n ${PER_POP} samples_by_pop/samples.yri.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.gbr.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.chb.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.clm.nofam.txt; } > samples.admix_clm.gbr.i35.txt

# ADMIX_CLM: YRI_IBS_CHB_CLM
{ head -n ${PER_POP} samples_by_pop/samples.yri.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.ibs.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.chb.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.clm.nofam.txt; } > samples.admix_clm.ibs.i35.txt

# ADMIX_PEL: YRI_CEU_CHB_PEL
{ head -n ${PER_POP} samples_by_pop/samples.yri.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.ceu.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.chb.nofam.txt &&
  head -n ${PER_POP} samples_by_pop/samples.pel.nofam.txt; } > samples.admix_pel.i35.txt
